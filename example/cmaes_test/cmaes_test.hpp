//
// Created by yyk-22 on 24-1-1.
//

#ifndef CMAES_TEST_CMAES_TEST_HPP
#define CMAES_TEST_CMAES_TEST_HPP
#include "cmaes_interface.h"

#include "SEvoBench/single_algorithm.hpp"
#include <cstring>
template <typename T, int Dim, int Pop_Size> class cmaes_wrapper {
  std::array<T, Dim> initial_solution;
  std::array<T, Dim> initial_sigma;
  T *const *X, *y;
  cmaes_t cma;
  T l, r;

public:
  cmaes_wrapper(T _l, T _r, T sigma_rate) : l(_l), r(_r) {
    initial_solution.fill(_l + (_r - _l) / T(2));
    initial_sigma.fill((_r - _l) * sigma_rate);
    std::random_device rd{};
    cmaes_init_para(&cma, Dim, initial_solution.data(), initial_sigma.data(),
                    rd(), Pop_Size, "no");
    cma.sp.filename = strdup("no");
    y = cmaes_init_final(&cma);
  }
  void sample_population() { X = cmaes_SamplePopulation(&cma); }
  auto get_eval() { return cmaes_Get(&cma, "eval"); }
  auto get_fbestever() { return cmaes_Get(&cma, "fbestever"); }
  auto box_handle(int i) {
    constexpr int sample_max = 20;
    for (int ite = 0; std::any_of(X[i], X[i] + Dim,
                                  [this](auto x) { return x > r || x < l; });
         ite++) {
      cmaes_ReSampleSingle(&cma, i);
      if (ite > sample_max) {
        for (int k = 0; k < Dim; k++)
          X[i][k] = (l + r) / T(2);
        return;
      }
    }
  }
  auto update_distribution() { cmaes_UpdateDistribution(&cma, y); }
  template <typename F> void set(F &&f, int i) { y[i] = f(X[i]); }
  auto get_xbestever() {
    std::array<T, Dim> result;
    std::copy_n(cmaes_GetPtr(&cma, "xbestever"), Dim, result.data());
    return result;
  }
  ~cmaes_wrapper() { cmaes_exit(&cma); }
};

template <typename T> struct cmaes_parameter {
  T sigma_rate = T(1) / T(6);
};
template <int Dim, int Pop_Size, int Max, int Memory_Flag, typename F,
          typename T, typename Parameter_Type = cmaes_parameter<T>>
inline auto cmaes(F &&f, T l, T r, const Parameter_Type &pt = Parameter_Type()) {
  sevobench::tool::curve_t<Memory_Flag, T, Max> cmaes_curve{};
  cmaes_wrapper<T, Dim, Pop_Size> cma_wrap(l, r, pt.sigma_rate);
  for (int _ = 0; (Memory_Flag ? _ : cma_wrap.get_eval()) < Max;
       _ += (Memory_Flag ? 1 : 0)) {
    cma_wrap.sample_population();
    for (size_t i = 0; i < Pop_Size; ++i) {
      cma_wrap.box_handle(i);
      cma_wrap.set(f, i);
    }
    cma_wrap.update_distribution();
    if constexpr (Memory_Flag) {
      cmaes_curve[2 * _] = cma_wrap.get_eval();
      cmaes_curve[2 * _ + 1] = cma_wrap.get_fbestever();
    }
  }
  if constexpr (Memory_Flag)
    return std::make_tuple(cma_wrap.get_xbestever(), cma_wrap.get_fbestever(),
                           std::move(cmaes_curve));
  if constexpr (!Memory_Flag)
    return std::make_tuple(cma_wrap.get_xbestever(), cma_wrap.get_fbestever());
}
REGISTER_ALGORITHM(CMA_ES, cmaes, cmaes_parameter);
#endif // CMAES_TEST_CMAES_TEST_HPP
