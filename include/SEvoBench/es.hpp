#pragma once
#include "tool.hpp"

namespace sevobench {

template <std::floating_point T> struct es_parameter {
  T sigma = T(1);
};

template <typename P, typename T>
concept es_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.sigma)>>;
};

template <int Dim, int Pop_Size = 1, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename G, typename F,
          std::floating_point T, typename Parameter_Type = es_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           es_parameter_concept<Parameter_Type, T> &&
           algorithm_vector_concept<F, G, T>
inline auto es_optimize(G &&sol, F &&function, T l, T r,
                        const Parameter_Type &pt = Parameter_Type()) noexcept {
  T alpha = std::pow(T(2), T(1) / T(Dim));
  T beta = std::pow(T(2), T(-0.25) / T(Dim));
  T sigma = pt.sigma;
  std::random_device rd;
  tool::simple_rand sr1(rd()), sr2(rd());
  std::remove_reference_t<G> tmp;
  tool::curve_t<Memory_Flag, T, Max> es_convergence_curve;
  auto fit = function(sol.data());
  constexpr T k = T(1) / T(0x7fff);
  if constexpr (!Memory_Flag) {
    for (int i = 1; i < Max; i++) {
      for (int j = 0; j < Dim; j++)
        tmp[j] = std::clamp(
            sol[j] + tool::box_muller(T(0), sigma, sr1() * k, sr2() * k), l, r);
      auto y = function(tmp.data());
      if (y < fit) {
        sol = tmp;
        fit = y;
        sigma *= alpha;
      } else
        sigma *= beta;
    }
  } else {
    for (int i = 0; i < Max; i++) {
      for (int _ = 0; _ < Pop_Size; _++) {
        for (int j = 0; j < Dim; j++)
          tmp[j] = std::clamp(
              sol[j] + tool::box_muller(T(0), sigma, sr1() * k, sr2() * k), l,
              r);
        auto y = function(tmp.data());
        if (y < fit) {
          sol = tmp;
          fit = y;
          sigma *= alpha;
        } else
          sigma *= beta;
      }
      es_convergence_curve[2 * i] = (i + 1) * Pop_Size + 1;
      es_convergence_curve[2 * i + 1] = fit;
    }
  }
  if constexpr (Memory_Flag)
    return std::make_tuple(sol, fit, std::move(es_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(sol, fit);
}

template <int Dim, int Pop_Size = 1, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = es_parameter<T>>
inline auto es(F &&function, T l, T r,
               const Parameter_Type &pt = Parameter_Type()) noexcept {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<T> dis(l, r);
  std::array<T, Dim> sol;
  for (auto &x : sol)
    x = dis(gen);
  return es_optimize<Dim, Pop_Size, Max, Memory_Flag>(sol, function, l, r, pt);
}
} // namespace sevobench
