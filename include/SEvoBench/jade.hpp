#pragma once

#include "tool.hpp"

namespace sevobench {

template <std::floating_point T> struct jade_parameter {
  T c = 0.1;
  T p = 0.1;
};

template <typename P, typename T>
concept jade_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.c)>>;
  requires std::same_as<T, std::remove_cvref_t<decltype(p.p)>>;
};

template <int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename G, typename F,
          std::floating_point T, typename Parameter_Type = jade_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           jade_parameter_concept<Parameter_Type, T> &&
           algorithm_positions_concept<F, G, T> && (Pop_Size >= 3)
inline auto jade_optimize(
    G &&positions, F &&function, T l, T r,
    const Parameter_Type &pt = Parameter_Type()) noexcept {
  const T c = pt.c;
  const int P = static_cast<int>(Pop_Size * pt.p);
  constexpr T k = T(1) / T(0x7fff);
  T cr = 0.5;
  T f = 0.5;
  std::random_device rd;
  tool::simple_rand sr1(rd()), sr2(rd()), sr3(rd()), sr4(rd()), sr5(rd());
  std::vector<std::array<T, Dim>> tmp(Pop_Size);
  std::vector<std::array<T, Dim>> Archive(Pop_Size);
  int Archive_Size = 0;
  int replace_index = 0;
  tool::curve_t<Memory_Flag, T, Max> jade_convergence_curve;
  std::array<T, Pop_Size> fit;
  std::array<T, Pop_Size> tmp_fit;
  std::array<int, Pop_Size> index;
  std::array<T, Dim> random_c;
  std::iota(index.begin(), index.end(), 0);
  std::transform(positions.begin(), positions.begin() + Pop_Size, fit.begin(),
                 [&](auto &x) { return function(x.data()); });
  auto gen_f = [](auto x, auto &g) {
    if (!std::isfinite(x))
      return T(1);
    auto z = tool::cauchy_dis(x, T(0.1), g() * k);
    while (z <= 0)
      z = tool::cauchy_dis(x, T(0.1), g() * k);
    return z < 1 ? z : 1;
  };
  for (int _ = (Memory_Flag ? 0 : Pop_Size); _ < Max;
       _ += (Memory_Flag ? 1 : Pop_Size)) {
    std::array<T, Pop_Size> f_archive;
    std::array<T, Pop_Size> cr_archive;
    auto f_first = f_archive.begin();
    auto cr_first = cr_archive.begin();
    std::nth_element(index.begin(), index.begin() + P, index.end(),
                     [&](auto x, auto y) { return fit[x] < fit[y]; });
    for (int i = 0; i < Pop_Size; i++) {
      auto cR = std::clamp(tool::box_muller(cr, T(0.1), sr1() * k, sr2() * k),
                           T(0), T(1));
      auto cF = gen_f(f, sr3);
      int r0, r1, r2;
      r0 = index[P == 0 ? 0 : sr4() % P];
      do {
        r1 = sr4() % Pop_Size;
      } while (r1 == i);
      do {
        r2 = sr4() % (Pop_Size + Archive_Size);
      } while (r2 == i || r2 == r1);
      auto &vr2 = r2 >= Pop_Size ? Archive[r2 - Pop_Size] : positions[r2];
      auto j_rand = sr4() % Dim;
      std::generate_n(random_c.begin(), Dim, [&] { return sr5() * k; });

      for (int j = 0; j < Dim; j++) {
        auto x = (j == j_rand || random_c[j] < cR)
                     ? positions[i][j] +
                           cF * (positions[r0][j] - positions[i][j]) +
                           cF * (positions[r1][j] - vr2[j])
                     : positions[i][j];
        if (x < l)
          x = (positions[i][j] + l) * T(0.5);
        if (x > r)
          x = (positions[i][j] + r) * T(0.5);
        tmp[i][j] = x;
      }
      tmp_fit[i] = function(tmp[i].data());
      if (tmp_fit[i] < fit[i]) {
        *f_first++ = cF;
        *cr_first++ = cR;
      }
    }
    for (int i = 0; i < Pop_Size; i++)
      if (tmp_fit[i] < fit[i]) {
        fit[i] = tmp_fit[i];
        positions[i] = tmp[i];
        Archive[replace_index++] = positions[i];
        replace_index %= Pop_Size;
        Archive_Size = Archive_Size == Pop_Size ? Pop_Size : Archive_Size + 1;
      }
    cr = (1 - c) * cr +
         (cr_first == cr_archive.begin()
              ? 0
              : c * std::accumulate(cr_archive.begin(), cr_first, T(0)) /
                    (cr_first - cr_archive.begin()));
    f = (1 - c) * f +
        (f_first == f_archive.begin()
             ? 0
             : c *
                   std::inner_product(f_archive.begin(), f_first,
                                      f_archive.begin(), T(0)) /
                   std::accumulate(f_archive.begin(), f_first, T(0)));
    if constexpr (Memory_Flag) {
      jade_convergence_curve[2 * _] = (_ + 2) * Pop_Size;
      jade_convergence_curve[2 * _ + 1] =
          *std::min_element(fit.begin(), fit.end());
    }
  }
  auto i = std::min_element(fit.begin(), fit.end()) - fit.begin();
  if constexpr (Memory_Flag)
    return std::make_tuple(positions[i], fit[i],
                           std::move(jade_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(positions[i], fit[i]);
}

template <int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = jade_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           jade_parameter_concept<Parameter_Type, T> && (Pop_Size >= 3)
inline auto jade(F &&function, T l, T r,
                 const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return jade_optimize<Dim, Pop_Size, Max, Memory_Flag>(pop, function, l, r,
                                                        pt);
}

} // namespace sevobench
