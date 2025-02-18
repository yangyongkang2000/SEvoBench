#pragma once

#include "../common/tool.hpp"
#include "problem.hpp"
namespace sevobench::problem {

template <template <int, int, typename> class Problem, int Dim,
          std::floating_point T>
class suite {
private:
  std::vector<std::unique_ptr<single_problem<T>>> problems;
  std::vector<int> prob_ids;
  int ins_count = 0;

public:
  using value_type = T;
  static constexpr auto dim() { return Dim; }
  auto begin() const { return problems.begin(); }
  auto end() const { return problems.end(); }
  auto size() const { return static_cast<int>(problems.size()); }
  auto instance_count() const { return ins_count; }
  auto problem_index() const { return prob_ids; }
  template <typename Arg>
  suite(std::vector<int> arg1, Arg arg2)
      : problems(detail::create_problems<Problem, Dim, T>(arg1, arg2)),
        prob_ids(detail::filter_problem_indexs<Problem, Dim, T>(arg1)),
        ins_count([&] {
          if constexpr (std::is_same_v<Arg, int>) {
            return std::max(arg2, 1);
          } else
            return 1;
        }()) {}
};

namespace detail {
struct suite_args {
  std::vector<int> arg1;
  std::string dir_name;
  int ins_count;
};
} // namespace detail

template <template <int, int, typename> class Problem, int Dim = 0,
          std::floating_point T = float, bool B2 = false, bool B3 = true,
          bool B4 = false, bool B5 = false, bool B6 = false>
class suite_builder {
  detail::suite_args args;

public:
  suite_builder(const detail::suite_args &_args) : args(_args) {}
  suite_builder(detail::suite_args &&_args) : args(std::move(_args)) {}
  suite_builder() = default;
  template <int N> auto dim() {
    return suite_builder<Problem, N, T, true, B3, B4, B5, B6>(std::move(args));
  }
  template <std::floating_point T1> auto type() {
    return suite_builder<Problem, Dim, T1, B2, true, B4, B5, B6>(
        std::move(args));
  }
  auto dir(std::string dir) {
    args.dir_name = std::move(dir);
    return suite_builder<Problem, Dim, T, B2, B3, B4, true, false>(
        std::move(args));
  }
  auto instance_count(int count) {
    args.ins_count = count;
    return suite_builder<Problem, Dim, T, B2, B3, B4, false, true>(
        std::move(args));
  }
  auto problem_index(std::vector<int> prob_ids) {
    args.arg1 = std::move(prob_ids);
    return suite_builder<Problem, Dim, T, B2, B3, true, B5, B6>(
        std::move(args));
  }
  auto build() {
    static_assert(B2, "NO SET DIM");
    static_assert(B3, "NO SET TYPE");
    static_assert(B4, "NO SET PROBLEM INDEXS");
    static_assert(B5 || B6, "NO SET INSTANCE COUNT/DIR");
    if constexpr (B5) {
      return suite<Problem, Dim, T>(std::move(args.arg1),
                                    std::move(args.dir_name));
    } else {
      return suite<Problem, Dim, T>(std::move(args.arg1), args.ins_count);
    }
  }
};
} // namespace sevobench::problem
