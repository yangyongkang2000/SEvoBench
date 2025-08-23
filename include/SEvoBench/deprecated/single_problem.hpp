#pragma once
#include "cec2020.hpp"
#include "cec2022.hpp"
#include "classic_problem.hpp"
#include "tool.hpp"

#define SEvoBench_REGISTER_PROBLEMS(problem_name, problem, problem_size)       \
  namespace sevobench {                                                        \
  static constexpr auto problem_name =                                         \
      tool::alg_hash<std::uint64_t>("sevobench::" #problem);                   \
  template <int Prob_Index, int Dim, std::floating_point T>                    \
    requires(problem_size >= 1 && Prob_Index >= 1 &&                           \
             Prob_Index <= problem_size && Dim >= 1)                           \
  struct single_problem<problem_name, Prob_Index, Dim, T> {                    \
    static constexpr char name[] = #problem_name;                              \
    static constexpr int size = problem_size;                                  \
    static constexpr int index = Prob_Index;                                   \
    static constexpr int dim = Dim;                                            \
    static constexpr T L = problem<Prob_Index, Dim, T>::L;                     \
    static constexpr T U = problem<Prob_Index, Dim, T>::U;                     \
    T operator()(T *x) noexcept { return problem<Prob_Index, Dim, T>()(x); }   \
    T operator()(const T *x, int) noexcept { return (*this)(x); }              \
    using value_type = T;                                                      \
  };                                                                           \
  }

#define SEvoBench_INIT_PROBLEMS(problem, F)                                    \
  namespace sevobench {                                                        \
  template <int Dim, std::floating_point T>                                    \
    requires(Dim >= 1)                                                         \
  struct init_single_problem<                                                  \
      tool::alg_hash<std::uint64_t>("sevobench::" #problem), Dim, T> {         \
    inline static bool is_init = false;                                        \
    init_single_problem() {                                                    \
      F<T, Dim>();                                                             \
      is_init = true;                                                          \
    }                                                                          \
  };                                                                           \
  }

namespace sevobench {

template <std::uint64_t Prob_HashName, int Prob_Index, int Dim,
          std::floating_point T>
struct single_problem;

template <std::uint64_t Prob_HashName, int Dim, std::floating_point T>
struct init_single_problem {
  inline static bool is_init = true;
};

constexpr auto ShiftFunc =
    tool::alg_hash<std::uint64_t>("sevobench::ShiftFunc");
constexpr auto RotateShiftFunc =
    tool::alg_hash<std::uint64_t>("sevobench::RotateShiftFunc");
constexpr auto AsyShiftFunc =
    tool::alg_hash<std::uint64_t>("sevobench::AsyShiftFunc");
constexpr auto AsyRotateShiftFunc =
    tool::alg_hash<std::uint64_t>("sevobench::AsyRotateShiftFunc");
constexpr auto CEC2020 =
    tool::alg_hash<std::uint64_t>("sevobench::ieee_cec_set::cec2020_problem");
constexpr auto CEC2022 =
    tool::alg_hash<std::uint64_t>("sevobench::ieee_cec_set::cec2022_problem");

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct init_single_problem<ShiftFunc, Dim, T> {
  inline static bool is_init = false;
  init_single_problem() { yao_func::init_shift_rotate<false, Dim, T>(); }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct init_single_problem<AsyShiftFunc, Dim, T> {
  inline static bool is_init = false;
  init_single_problem() { yao_func::init_shift_rotate<false, Dim, T>(); }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct init_single_problem<RotateShiftFunc, Dim, T> {
  inline static bool is_init = false;
  init_single_problem() { yao_func::init_shift_rotate<true, Dim, T>(); }
};
template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct init_single_problem<AsyRotateShiftFunc, Dim, T> {
  inline static bool is_init = false;
  init_single_problem() { yao_func::init_shift_rotate<true, Dim, T>(); }
};

template <int Prob_Index, int Dim, std::floating_point T>
  requires(Prob_Index >= 1 && Prob_Index <= 13 && Dim >= 1)
struct single_problem<ShiftFunc, Prob_Index, Dim, T> {
  static constexpr char name[] = "ShiftFunc";
  static constexpr int size = 13;
  static constexpr int index = Prob_Index;
  static constexpr int dim = Dim;
  static constexpr T L =
      yao_func::classic_problem<Prob_Index, Dim, T, false, false>::L;
  static constexpr T U =
      yao_func::classic_problem<Prob_Index, Dim, T, false, false>::U;
  T operator()(const T *x) noexcept {
    return yao_func::classic_problem<Prob_Index, Dim, T, false, false>()(
        x, yao_func::detail::shift_vector_data<T, Dim>[Prob_Index - 1],
        yao_func::detail::rotate_matrix_data<T, Dim>[Prob_Index - 1]);
  }
  T operator()(const T *x, int) noexcept { return (*this)(x); }
  using value_type = T;
};

template <int Prob_Index, int Dim, std::floating_point T>
  requires(Prob_Index >= 1 && Prob_Index <= 13 && Dim >= 1)
struct single_problem<AsyShiftFunc, Prob_Index, Dim, T> {
  static constexpr char name[] = "AsyShiftFunc";
  static constexpr int size = 13;
  static constexpr int index = Prob_Index;
  static constexpr int dim = Dim;
  static constexpr T L =
      yao_func::classic_problem<Prob_Index, Dim, T, false, true>::L;
  static constexpr T U =
      yao_func::classic_problem<Prob_Index, Dim, T, false, true>::U;
  T operator()(const T *x) noexcept {
    return yao_func::classic_problem<Prob_Index, Dim, T, false, true>()(
        x, yao_func::detail::shift_vector_data<T, Dim>[Prob_Index - 1],
        yao_func::detail::rotate_matrix_data<T, Dim>[Prob_Index - 1]);
  }
  T operator()(const T *x, int) noexcept { return (*this)(x); }
  using value_type = T;
};

template <int Prob_Index, int Dim, std::floating_point T>
  requires(Prob_Index >= 1 && Prob_Index <= 13 && Dim >= 1)
struct single_problem<RotateShiftFunc, Prob_Index, Dim, T> {
  static constexpr char name[] = "RotateShiftFunc";
  static constexpr int size = 13;
  static constexpr int index = Prob_Index;
  static constexpr int dim = Dim;
  static constexpr T L =
      yao_func::classic_problem<Prob_Index, Dim, T, true, false>::L;
  static constexpr T U =
      yao_func::classic_problem<Prob_Index, Dim, T, true, false>::U;
  T operator()(const T *x) noexcept {
    return yao_func::classic_problem<Prob_Index, Dim, T, true, false>()(
        x, yao_func::detail::shift_vector_data<T, Dim>[Prob_Index - 1],
        yao_func::detail::rotate_matrix_data<T, Dim>[Prob_Index - 1]);
  }
  T operator()(const T *x, int) noexcept { return (*this)(x); }
  using value_type = T;
};

template <int Prob_Index, int Dim, std::floating_point T>
  requires(Prob_Index >= 1 && Prob_Index <= 13 && Dim >= 1)
struct single_problem<AsyRotateShiftFunc, Prob_Index, Dim, T> {
  static constexpr char name[] = "AsyRotateShiftFunc";
  static constexpr int size = 13;
  static constexpr int index = Prob_Index;
  static constexpr int dim = Dim;
  static constexpr T L =
      yao_func::classic_problem<Prob_Index, Dim, T, true, true>::L;
  static constexpr T U =
      yao_func::classic_problem<Prob_Index, Dim, T, true, true>::U;
  T operator()(const T *x) noexcept {
    return yao_func::classic_problem<Prob_Index, Dim, T, true, true>()(
        x, yao_func::detail::shift_vector_data<T, Dim>[Prob_Index - 1],
        yao_func::detail::rotate_matrix_data<T, Dim>[Prob_Index - 1]);
  }
  T operator()(const T *x, int) noexcept { return (*this)(x); }
  using value_type = T;
};
template <int Prob_Index, int Dim, std::floating_point T>
  requires(Prob_Index >= 1 && Prob_Index <= 10 && Dim >= 1)
struct single_problem<CEC2020, Prob_Index, Dim, T> {
  static constexpr char name[] = "CEC2020";
  static constexpr int size = 10;
  static constexpr int index = Prob_Index;
  static constexpr int dim = Dim;
  static constexpr T L = ieee_cec_set::cec2020_problem<Prob_Index, Dim, T>::L;
  static constexpr T U = ieee_cec_set::cec2020_problem<Prob_Index, Dim, T>::U;
  T operator()(const T *x) noexcept {
    if constexpr (Prob_Index < 8) {
      return ieee_cec_set::cec2020_problem<Prob_Index, Dim, T>()(
          x,
          ieee_cec_set::cec2020_data::shift_vector_data<T, Dim>[Prob_Index - 1],
          ieee_cec_set::cec2020_data::rotate_matrix_data<T, Dim>[Prob_Index -
                                                                 1]);
    }
    if constexpr (Prob_Index > 7) {
      return ieee_cec_set::cec2020_problem<Prob_Index, Dim, T>()(x);
    }
  }
  T operator()(const T *x, int) noexcept { return (*this)(x); }
  using value_type = T;
};

template <int Prob_Index, int Dim, std::floating_point T>
  requires(Prob_Index >= 1 && Prob_Index <= 12 && Dim >= 1)
struct single_problem<CEC2022, Prob_Index, Dim, T> {
  static constexpr char name[] = "CEC2022";
  static constexpr int size = 12;
  static constexpr int index = Prob_Index;
  static constexpr int dim = Dim;
  static constexpr T L = ieee_cec_set::cec2022_problem<Prob_Index, Dim, T>::L;
  static constexpr T U = ieee_cec_set::cec2022_problem<Prob_Index, Dim, T>::U;
  T operator()(const T *x) noexcept {
    if constexpr (Prob_Index < 9) {
      return ieee_cec_set::cec2022_problem<Prob_Index, Dim, T>()(
          x,
          ieee_cec_set::cec2022_data::shift_vector_data<T, Dim>[Prob_Index - 1],
          ieee_cec_set::cec2022_data::rotate_matrix_data<T, Dim>[Prob_Index -
                                                                 1]);
    }
    if constexpr (Prob_Index > 8) {
      return ieee_cec_set::cec2022_problem<Prob_Index, Dim, T>()(x);
    }
  }
  T operator()(const T *x, int) noexcept { return (*this)(x); }
  using value_type = T;
};

} // namespace sevobench

SEvoBench_INIT_PROBLEMS(ieee_cec_set::cec2020_problem,
                        ieee_cec_set::init_cec2020_problem)
    SEvoBench_INIT_PROBLEMS(ieee_cec_set::cec2022_problem,
                            ieee_cec_set::init_cec2022_problem)

#define REGISTER_PROBLEMS(problem_name, problem, problem_size)                 \
  static constexpr auto problem_name =                                         \
      sevobench::tool::alg_hash<std::uint64_t>(#problem);                      \
  template <int Prob_Index, int Dim, std::floating_point T>                    \
    requires(problem_size >= 1 && Prob_Index >= 1 &&                           \
             Prob_Index <= problem_size && Dim >= 1)                           \
  struct sevobench::single_problem<problem_name, Prob_Index, Dim, T> {         \
    static constexpr char name[] = #problem_name;                              \
    static constexpr int size = problem_size;                                  \
    static constexpr int index = Prob_Index;                                   \
    static constexpr int dim = Dim;                                            \
    static constexpr T L = problem<Prob_Index, Dim, T>::L;                     \
    static constexpr T U = problem<Prob_Index, Dim, T>::U;                     \
    T operator()(T *x) noexcept { return problem<Prob_Index, Dim, T>()(x); }   \
    T operator()(const T *x, int) noexcept { return (*this)(x); }              \
    using value_type = T;                                                      \
  };

#define INIT_PROBLEMS(problem, F)                                              \
  template <int Dim, std::floating_point T>                                    \
    requires(Dim >= 1)                                                         \
  struct sevobench::init_single_problem<                                       \
      sevobench::tool::alg_hash<std::uint64_t>(#problem), Dim, T> {            \
    inline static bool is_init = false;                                        \
    init_single_problem() {                                                    \
      F<T, Dim>();                                                             \
      is_init = true;                                                          \
    }                                                                          \
  };
