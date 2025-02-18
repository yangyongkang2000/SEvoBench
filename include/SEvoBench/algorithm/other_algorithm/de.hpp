#pragma once
#include "../../common/common_concept.hpp"
#include "../../common/tool.hpp"

namespace sevobench::other_algorithm {

template <std::floating_point T> struct de_parameter {
  T f = 0.5;
  T cr = 0.9;
};

template <typename P, typename T>
concept de_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.f)>>;
  requires std::same_as<T, std::remove_cvref_t<decltype(p.cr)>>;
};

template <int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename G, typename F,
          std::floating_point T, typename Parameter_Type = de_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           de_parameter_concept<Parameter_Type, T> &&
           algorithm_positions_concept<F, G, T> && (Pop_Size >= 4)
inline auto de_optimize(G &&positions, F &&function, T l, T r,
                        const Parameter_Type &pt = Parameter_Type()) noexcept {
  const T f = pt.f;
  const int c = static_cast<int>(pt.cr * 0x7fff);
  std::vector<std::array<T, Dim>> tmp(Pop_Size);
  tool::curve_t<Memory_Flag, T, Max> de_convergence_curve;
  std::array<T, Pop_Size> fit;
  std::array<T, Pop_Size> tmp_fit;
  std::array<int, Dim> random_c;
  std::random_device rd;
  tool::simple_rand sr1(rd()), sr2(rd());
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<T> dis(l, r);
  std::transform(positions.begin(), positions.begin() + Pop_Size, fit.begin(),
                 [&](auto &x) { return function(x.data()); });
  for (int _ = (Memory_Flag ? 0 : Pop_Size); _ < Max;
       _ += (Memory_Flag ? 1 : Pop_Size)) {
    for (int i = 0; i < Pop_Size; i++) {
      int r1, r2, r3;
      std::generate_n(random_c.begin(), Dim, [&] { return sr1(); });
      do {
        r1 = sr2() % Pop_Size;
      } while (r1 == i);
      do {
        r2 = sr2() % Pop_Size;
      } while (r2 == i || r2 == r1);
      do {
        r3 = sr2() % Pop_Size;
      } while (r3 == i || r3 == r1 || r3 == r2);
      for (int j = 0; j < Dim; j++) {
        tmp[i][j] = random_c[j] < c
                        ? std::clamp(positions[r1][j] + f * (positions[r2][j] -
                                                             positions[r3][j]),
                                     l, r)
                        : positions[i][j];
      }
      tmp_fit[i] = function(tmp[i].data());
    }
    for (int i = 0; i < Pop_Size; i++)
      if (tmp_fit[i] < fit[i]) {
        fit[i] = tmp_fit[i];
        std::copy_n(tmp[i].data(), Dim, positions[i].data());
      }
    if constexpr (Memory_Flag) {
      de_convergence_curve[2 * _] = (_ + 2) * Pop_Size;
      de_convergence_curve[2 * _ + 1] =
          *std::min_element(fit.begin(), fit.end());
    }
  }
  auto index = std::min_element(fit.begin(), fit.end()) - fit.begin();
  if constexpr (Memory_Flag)
    return std::make_tuple(positions[index], fit[index],
                           std::move(de_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(positions[index], fit[index]);
}

template <int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = de_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           de_parameter_concept<Parameter_Type, T> && (Pop_Size >= 4)
inline auto de(F &&function, T l, T r,
               const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return de_optimize<Dim, Pop_Size, Max, Memory_Flag>(pop, function, l, r, pt);
}

} // namespace sevobench::other_algorithm
