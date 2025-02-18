#pragma once

#include "../problem/suite.hpp"
#include "../utility/parallel_task.hpp"

namespace sevobench::experiment {

template <std::floating_point T> struct suite_observer {
  virtual void log(const problem::problem_state<T> &,
                   const problem::problem_info<T> &) = 0;
  virtual ~suite_observer() = default;
};

namespace detail {
template <std::floating_point T, typename S> struct suite_problem {
private:
  problem::single_problem<T> *p;
  S &o;
  mutable int run_id = 0;
  mutable int evals = 0;

public:
  suite_problem(problem::single_problem<T> *_p, S &_o, int _run_id)
      : p(_p), o(_o), run_id(_run_id) {}
  auto operator()(std::span<const T> x) const noexcept {
    auto value = (*p)(x);
    o.log(problem::problem_state<T>{.evaluations = ++evals,
                                    .run_id = run_id,
                                    .current_value = value,
                                    .current_x = x},
          p->problem_information());
    return value;
  }
  auto lower_bound() const noexcept { return p->lower_bound(); }
  auto upper_bound() const noexcept { return p->upper_bound(); }
};

} // namespace detail

template <bool parallel = true>
inline auto evo_bench(auto &&alg, auto &&su, auto &&obs,
                      int independent_runs) noexcept {
  [[maybe_unused]] std::conditional_t<parallel, parallel_task, void *> pt{};
  [[maybe_unused]] std::conditional_t<parallel, std::vector<std::future<void>>,
                                      void *> v{};
  if constexpr (parallel) {
    v.reserve(su.size() * independent_runs);
  }
  for (const auto &p : su) {
    constexpr auto b1 =
        requires { alg(*p, p->lower_bound(), p->upper_bound()); };
    constexpr auto b2 = requires { alg(*p); };
    static_assert(b1 || b2, "ALGORITHM IS INVALID!");
    for (int i = 0; i < independent_runs; i++) {
      auto temp_p = detail::suite_problem(p.get(), obs, i + 1);
      if constexpr (b1) {
        if constexpr (parallel) {
          v.emplace_back(pt.submit([=] {
            alg(temp_p, temp_p.upper_bound(), temp_p.upper_bound());
          }));
        } else {
          alg(temp_p, temp_p.lower_bound(), temp_p.upper_bound());
        }
      } else {
        if constexpr (parallel) {
          v.emplace_back(pt.submit([=] { alg(temp_p); }));
        } else {
          alg(temp_p);
        }
      }
    }
  }
  if constexpr (parallel) {
    for (auto &_ : v)
      _.get();
  }
}
} // namespace sevobench::experiment
