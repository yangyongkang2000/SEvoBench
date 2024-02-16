#pragma once

#include "tool.hpp"

namespace sevobench {

template <std::floating_point T> struct cso_parameter {
  T fai = T(0.1);
};

template <typename P, typename T>
concept cso_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.fai)>>;
};

template <int Dim, int Pop_Size = 250, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename G, typename F,
          std::floating_point T, typename Parameter_Type = cso_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           cso_parameter_concept<Parameter_Type, T> &&
           algorithm_positions_concept<F, G, T> && (Pop_Size >= 2)
inline auto cso_optimize(G &&positions, F &&function, T l, T r,
                         const Parameter_Type &pt = Parameter_Type()) noexcept {
  const T fai = pt.fai;
  constexpr T k = T(1) / T(0x7FFF);
  const T k1 = fai * k;
  constexpr T q = T(1) / T(Pop_Size);
  std::vector<std::array<T, Dim>> vec(Pop_Size);
  tool::curve_t<Memory_Flag, T, Max> cso_convergence_curve;
  std::array<T, Pop_Size> fit;
  std::array<T, Dim> r1;
  std::array<T, Dim> r2;
  std::array<T, Dim> r3;
  std::array<int, Pop_Size> index;
  std::random_device rd;
  std::default_random_engine gen(rd()), seed_gen(rd());
  std::uniform_int_distribution<unsigned int> seed(0, 0x7fff);
  tool::simple_rand sr1(rd());
  tool::simple_rand sr2(rd());
  tool::simple_rand sr3(rd());
  std::transform(positions.begin(), positions.begin() + Pop_Size, fit.begin(),
                 [&](auto &x) { return function(x.data()); });
  std::iota(index.begin(), index.end(), 0);
  for (int _ = (Memory_Flag ? 0 : Pop_Size); _ < Max;
       _ += (Memory_Flag ? 1 : Pop_Size / 2)) {
    tool::kunth_shuffle(index.begin(), index.end(), seed(seed_gen));
    std::array<T, Dim> mean{};
    for (int i = 0; i < Pop_Size; i++)
      for (int j = 0; j < Dim; j++)
        mean[j] += q * positions[i][j];
    for (int i = 0; i < Pop_Size / 2; i++) {
      auto lo = index[i];
      auto w = index[i + Pop_Size / 2];
      if (fit[lo] < fit[w])
        std::swap(lo, w);
      for (int j = 0; j < Dim; j++) {
        r1[j] = sr1() * k;
        r2[j] = sr2() * k;
        r3[j] = sr3() * k1;
      }
      for (int j = 0; j < Dim; j++)
        vec[lo][j] = r1[j] * vec[lo][j] +
                     r2[j] * (positions[w][j] - positions[lo][j]) +
                     r3[j] * (mean[j] - positions[lo][j]);
      for (int j = 0; j < Dim; j++)
        positions[lo][j] = std::clamp(positions[lo][j] + vec[lo][j], l, r);
      fit[lo] = function(positions[lo].data());
    }
    if constexpr (Memory_Flag) {
      cso_convergence_curve[2 * _] = Pop_Size / 2 * (_ + 1) + Pop_Size;
      cso_convergence_curve[2 * _ + 1] =
          *std::min_element(fit.begin(), fit.end());
    }
  }
  auto i = std::min_element(fit.begin(), fit.end()) - fit.begin();
  if constexpr (Memory_Flag)
    return std::make_tuple(positions[i], fit[i],
                           std::move(cso_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(positions[i], fit[i]);
}

template <int Dim, int Pop_Size = 250, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = cso_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           cso_parameter_concept<Parameter_Type, T> && (Pop_Size >= 2)
inline auto cso(F &&function, T l, T r,
                const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return cso_optimize<Dim, Pop_Size, Max, Memory_Flag>(pop, function, l, r, pt);
}

} // namespace sevobench
