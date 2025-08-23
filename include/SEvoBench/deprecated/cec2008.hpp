#pragma once
#include "cec_basical_func_simd.hpp"
#include "single_problem.hpp"
namespace sevobench::ieee_cec_set::cec2008_data {
template <std::floating_point T, int Dim>
  requires simd_type_detail::simd_dim<T, Dim>
std::vector<std::array<T, Dim>> shift_vector_data(6);
}
namespace sevobench::ieee_cec_set {

template <std::floating_point T, int Dim>
  requires simd_type_detail::simd_dim<T, Dim>
inline void init_cec2008_problem() noexcept {
  std::random_device rd{};
  std::default_random_engine gen{rd()};
  static constexpr T L[] = {-100, -100, -101, -5, -600, -32};
  static constexpr T U[] = {100, 100, 99, 5, 600, 32};
  for (int i = 0; i < 6; i++) {
    std::uniform_real_distribution<T> dis(L[i], U[i]);
    std::generate_n(cec2008_data::shift_vector_data<T, Dim>[i].begin(), Dim,
                    [&] { return dis(gen); });
  }
}

template <int Prob_Index, int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim>
struct cec2008_problem {
  static constexpr T L =
      Prob_Index <= 3 ? -100
                      : (Prob_Index == 4 ? -5 : (Prob_Index == 5 ? -600 : -32));
  static constexpr T U =
      Prob_Index <= 3 ? 100
                      : (Prob_Index == 4 ? 5 : (Prob_Index == 5 ? 600 : 32));
  T operator()(const T *x) {
    std::array<T, Dim> tmp;
    CEC_Base_Func_Simd::shift_func<Dim>(
        x, tmp.data(),
        cec2008_data::shift_vector_data<T, Dim>[Prob_Index - 1].data());
    if constexpr (Prob_Index == 1) {
      return CEC_Base_Func_Simd::sphere<Dim>(tmp.data());
    }
    if constexpr (Prob_Index == 2) {
      return CEC_Base_Func_Simd::schwefel<Dim>(tmp.data());
    }
    if constexpr (Prob_Index == 3) {
      return CEC_Base_Func_Simd::rosenbrock<Dim>(tmp.data());
    }
    if constexpr (Prob_Index == 4) {
      return CEC_Base_Func_Simd::rastrigin<Dim>(tmp.data());
    }
    if constexpr (Prob_Index == 5) {
      return CEC_Base_Func_Simd::griewank<Dim>(tmp.data());
    }
    if constexpr (Prob_Index == 6) {
      return CEC_Base_Func_Simd::ackley<Dim>(tmp.data());
    }
  }
};
} // namespace sevobench::ieee_cec_set

SEvoBench_REGISTER_PROBLEMS(CEC2008, ieee_cec_set::cec2008_problem, 6)
    SEvoBench_INIT_PROBLEMS(ieee_cec_set::cec2008_problem,
                            ieee_cec_set::init_cec2008_problem)
