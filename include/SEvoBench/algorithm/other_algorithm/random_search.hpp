#pragma once
#include "../../common/common_concept.hpp"
#include "../../common/tool.hpp"

namespace sevobench::other_algorithm {

template <std::floating_point T> struct random_search_parameter {};

template <int Dim, int Pop_Size = 1, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename G, typename F,
          std::floating_point T>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           algorithm_vector_concept<F, G, T>
inline auto random_search_optimize(G &&sol, F &&function, T l, T r) noexcept {
  std::random_device rd;
  tool::simple_rand sr(rd());
  tool::curve_t<Memory_Flag, T, Max> random_search_convergence_curve;
  auto fit = function(sol.data());
  constexpr T k = T(1) / T(0x7fff);
  if constexpr (!Memory_Flag) {
    for (int i = 1; i < Max; i++) {
      auto index = sr() % Dim;
      auto tmp = sol[index];
      sol[index] = sr() * k * (r - l) + l;
      auto y = function(sol.data());
      if (y < fit)
        fit = y;
      else
        sol[index] = tmp;
    }
  } else {
    for (int i = 0; i < Max; i++) {
      for (int j = 0; j < Pop_Size; j++) {
        auto index = sr() % Dim;
        auto tmp = sol[index];
        sol[index] = sr() * k * (r - l) + l;
        auto y = function(sol.data());
        if (y < fit)
          fit = y;
        else
          sol[index] = tmp;
      }
      random_search_convergence_curve[2 * i] = (i + 1) * Pop_Size + 1;
      random_search_convergence_curve[2 * i + 1] = fit;
    }
  }
  if constexpr (Memory_Flag)
    return std::make_tuple(sol, fit,
                           std::move(random_search_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(sol, fit);
}

template <int Dim, int Pop_Size = 1, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = random_search_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T>
inline auto random_search(
    F &&function, T l, T r,
    [[maybe_unused]] const Parameter_Type &_ = Parameter_Type()) noexcept {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<T> dis(l, r);
  std::array<T, Dim> sol;
  for (auto &x : sol)
    x = dis(gen);
  return random_search_optimize<Dim, Pop_Size, Max, Memory_Flag>(sol, function,
                                                                 l, r);
}

} // namespace sevobench::other_algorithm
