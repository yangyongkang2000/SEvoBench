#pragma once
#include "../common/population.hpp"
#include "../common/tool.hpp"
namespace sevobench::problem {

template <int Index, int Dim, std::floating_point T>
  requires(Index > 0) && (Dim > 0)
struct problem_common {
  using value_type = T;
  constexpr static auto index() { return Index; }
  constexpr static auto dim() { return Dim; }
};

template <std::floating_point T> struct problem_info {
  int index;
  int instance;
  int dim;
  T lb;
  T ub;
  std::optional<solution<T>> optimum;
};

template <std::floating_point T> struct problem_state {
  int evaluations;
  int run_id;
  T current_value;
  std::span<const T> current_x;
};

template <int First_Index, int Last_Index>
  requires(First_Index > 0) && (First_Index <= Last_Index)
inline auto problem_range() noexcept {
  std::vector<int> v(Last_Index - First_Index + 1);
  std::iota(v.begin(), v.end(), First_Index);
  return v;
}

template <std::floating_point T> class single_problem {
  problem_info<T> problem_data_;

public:
  virtual T operator()(std::span<const T>) = 0;
  virtual ~single_problem() = default;
  auto index() { return problem_data_.index; }
  auto instance() { return problem_data_.instance; }
  auto dim() { return problem_data_.dim; }
  auto lower_bound() { return problem_data_.lb; }
  auto upper_bound() { return problem_data_.ub; }
  auto optimum() { return problem_data_.optimum; }
  auto optimum_num() {
    if (problem_data_.optimum)
      return std::optional<T>(problem_data_.optimum.fitness());
    else
      return std::optional<T>(std::nullopt);
  }
  const auto &problem_information() { return problem_data_; }
  single_problem() = default;
  single_problem(const single_problem &) = default;
  single_problem(single_problem &&) = default;
  single_problem &operator=(const single_problem &) = default;
  single_problem &operator=(single_problem &&) = default;
  single_problem(const problem_info<T> &_problem_data)
      : problem_data_(_problem_data) {}
  single_problem(problem_info<T> &&_problem_data)
      : problem_data_(std::move(_problem_data)) {}
};
namespace detail {

template <typename F>
concept single_problem_requires =
    requires(F &&f, std::span<const typename F::value_type> x) {
      requires std::is_floating_point_v<typename F::value_type>;
      {
        f.problem_information()
      } -> std::same_as<problem_info<typename F::value_type>>;
      F(std::move(f));
      { f(x) } -> std::same_as<typename F::value_type>;
    };

template <typename F>
  requires single_problem_requires<F>
class single_problem_wrapper : public single_problem<typename F::value_type> {
  F f;

public:
  using problem_type = F;
  using value_type = typename F::value_type;
  single_problem_wrapper(const single_problem_wrapper &) = default;
  single_problem_wrapper(single_problem_wrapper &&) = default;
  single_problem_wrapper &operator=(const single_problem_wrapper &) = default;
  single_problem_wrapper &operator=(single_problem_wrapper &&) = default;
  single_problem_wrapper(F &&_f)
      : single_problem<value_type>(_f.problem_information()), f(std::move(_f)) {
  }
  value_type operator()(std::span<const value_type> x) override { return f(x); }
};

template <template <int, int, typename> class Problem, int Dim,
          std::floating_point T, int I = 1, int S = 0>
inline constexpr auto get_problem_size() noexcept {
  if constexpr (single_problem_requires<Problem<I, Dim, T>>) {
    return get_problem_size<Problem, Dim, T, I + 1, S + 1>();
  } else {
    return S;
  }
}

template <template <int, int, typename> class Problem, int Dim,
          std::floating_point T, typename Arg, int I, int L>
  requires(L > 0) && requires(Arg arg) { Problem<I, Dim, T>(arg); }
inline auto generate_problem_factory(auto &v) noexcept {
  if constexpr (I <= L) {
    v.emplace_back([](Arg arg) {
      return std::unique_ptr<single_problem<T>>(
          std::make_unique<single_problem_wrapper<Problem<I, Dim, T>>>(
              Problem<I, Dim, T>(arg)));
    });
    if constexpr (I <= L - 1)
      generate_problem_factory<Problem, Dim, T, Arg, I + 1, L>(v);
  }
}

template <template <int, int, typename> class Problem, int Dim,
          std::floating_point T, typename Arg>
inline auto generate_problem_factory() noexcept {
  constexpr auto size = get_problem_size<Problem, Dim, T>();
  std::vector<std::unique_ptr<single_problem<T>> (*)(Arg)> v;
  v.reserve(size);
  generate_problem_factory<Problem, Dim, T, Arg, 1,
                           get_problem_size<Problem, Dim, T>()>(v);
  return v;
}

template <template <int, int, typename> class Problem, int Dim,
          std::floating_point T>
  requires(get_problem_size<Problem, Dim, T>() > 0)
inline auto filter_problem_indexs(std::vector<int> v) noexcept {
  constexpr auto size = get_problem_size<Problem, Dim, T>();
  if (!v.empty()) {
    for (auto &x : v)
      x = std::clamp(x, 1, size);
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
  }
  return v;
}

template <template <int, int, typename> class Problem, int Dim,
          std::floating_point T, typename Arg>
  requires requires { generate_problem_factory<Problem, Dim, T, Arg>(); }
inline auto create_problems(std::vector<int> arg1, Arg arg2) noexcept {
  std::vector<std::unique_ptr<single_problem<T>>> result;
  auto table = generate_problem_factory<Problem, Dim, T, Arg>();
  auto v = filter_problem_indexs<Problem, Dim, T>(arg1);
  if constexpr (std::is_same_v<Arg, int>) {
    arg2 = std::max(1, arg2);
  };
  if (!v.empty()) {
    if constexpr (!std::is_same_v<Arg, int>) {
      result.reserve(v.size());
    } else {
      result.reserve(arg2 * v.size());
    }
    for (auto index : v) {
      if constexpr (!std::is_same_v<Arg, int>) {
        result.emplace_back(table[index - 1](arg2));
      } else {
        for (int i = 1; i <= arg2; i++)
          result.emplace_back(table[index - 1](i));
      }
    }
  }
  return result;
}

} // namespace detail

} // namespace sevobench::problem
