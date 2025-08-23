#pragma once
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace sevobench {

namespace parallel_task_detail {
class task_function_wrapper {
  struct impl_base {
    virtual void call() = 0;
    virtual ~impl_base() {};
  };
  template <typename F> struct impl_type : impl_base {
    F f;
    void call() override { f(); }
    impl_type(F &&_f) : f(std::move(_f)) {}
  };
  std::unique_ptr<impl_base> impl;

public:
  void operator()() noexcept { impl->call(); }
  template <typename F>
  task_function_wrapper(F &&f) : impl(new impl_type<F>(std::move(f))) {}
  task_function_wrapper() = default;
  task_function_wrapper(task_function_wrapper &) = delete;
  task_function_wrapper(const task_function_wrapper &) = delete;
  task_function_wrapper &operator=(const task_function_wrapper &) = delete;
  task_function_wrapper(task_function_wrapper &&other)
      : impl(std::move(other.impl)) {}
  task_function_wrapper &operator=(task_function_wrapper &&other) {
    impl = std::move(other.impl);
    return *this;
  }
};
template <typename T> class threadsafe_queue {
  std::mutex mut;
  std::deque<T> data_queue;
  std::condition_variable data_cond;

public:
  void push(T new_value) noexcept {
    {
      std::lock_guard<std::mutex> lck(mut);
      data_queue.push_back(std::move(new_value));
    }
    data_cond.notify_one();
  }
  bool try_pop(T &value) noexcept {
    std::lock_guard<std::mutex> lck(mut);
    if (data_queue.empty())
      return false;
    value = std::move(data_queue.front());
    data_queue.pop_front();
    return true;
  }
};
} // namespace parallel_task_detail

class parallel_task {

  std::vector<std::thread> ts;
  parallel_task_detail::threadsafe_queue<
      parallel_task_detail::task_function_wrapper>
      tasks;
  bool done;
  void work_thread() noexcept {
    while (!done) {
      parallel_task_detail::task_function_wrapper task;
      if (tasks.try_pop(task)) {
        task();
      } else {
        std::this_thread::yield();
      }
    }
  }

public:
  parallel_task(unsigned int _sz) : done(false) {
    ts.reserve(_sz);
    for (unsigned int i = 0; i < _sz; i++)
      ts.emplace_back(&parallel_task::work_thread, this);
  }
  parallel_task() : parallel_task(std::thread::hardware_concurrency()) {}
  template <class Function, class... Args>
  [[nodiscard]] std::future<std::invoke_result_t<Function, Args...>>
  submit(Function &&f, Args &&...args) noexcept {
    auto g = [f = std::forward<Function>(f),
              ... args = std::forward<Args>(args)] { return f(args...); };
    auto pt(std::packaged_task<std::invoke_result_t<Function, Args...>()>(
        std::move(g)));
    auto res = pt.get_future();
    tasks.push(std::move(pt));
    return res;
  }
  ~parallel_task() noexcept {
    done = true;
    for (auto &_ : ts)
      if (_.joinable())
        _.join();
  }
};

} // namespace sevobench
