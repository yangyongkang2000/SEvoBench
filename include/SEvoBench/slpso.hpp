#pragma once
#include "tool.hpp"
namespace sevobench {

template <std::floating_point T> struct slpso_parameter {
  T beta = T(0.01);
  T alpha = T(0.5);
};

template <typename P, typename T>
concept slpso_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.beta)>>;
  requires std::same_as<T, std::remove_cvref_t<decltype(p.alpha)>>;
};

template <int Dim, int Pop_Size, int Max, bool Memory_Flag = false, typename G,
          typename F, std::floating_point T,
          typename Parameter_Type = slpso_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           slpso_parameter_concept<Parameter_Type, T> &&
           algorithm_positions_concept<F, G, T>
inline auto
slpso_optimize(G &&positions, F &&function, T l, T r,
               const Parameter_Type &pt = Parameter_Type()) noexcept {
  constexpr T q = T(1) / T(Pop_Size);
  const T alpha = pt.alpha;
  const T beta = pt.beta;
  const T eps = beta * T(Dim) / T(Pop_Size);
  constexpr T k = T(1) / T(0x7FFF);
  std::conditional_t<Dim <= Pop_Size, void *, std::array<T, Pop_Size>> P;
  std::conditional_t<Dim <= Pop_Size, void *, std::array<T, Pop_Size>> tmp_P;
  std::vector<std::array<T, Dim>> vec(Pop_Size);
  tool::curve_t<Memory_Flag, T, Max> slpso_convergence_curve;
  std::array<int, Pop_Size> index;
  std::array<int, Dim> win;
  std::array<T, Pop_Size> fit;
  std::array<T, Dim> r1;
  std::array<T, Dim> r2;
  std::array<T, Dim> r3;
  if constexpr (Dim > Pop_Size)
    for (int i = 0; i < Pop_Size - 1; i++)
      P[i] = std::pow(T(1) - q * T(i), alpha * std::log(std::ceil(T(Dim) * q)));
  std::iota(index.begin(), index.end(), 0);
  std::random_device rd;
  tool::simple_rand sr1(rd());
  tool::simple_rand sr2(rd());
  tool::simple_rand sr3(rd());
  tool::simple_rand sr4(rd());
  tool::simple_rand sr5(rd());
  std::transform(positions.begin(), positions.begin() + Pop_Size, fit.begin(),
                 [&](auto &x) { return function(x.data()); });
  T min_value = *std::min_element(fit.begin(), fit.end());
  int fes(Pop_Size);
  for (int _ = 0; (Memory_Flag ? _ : fes) < Max; _++) {
    std::array<T, Dim> mean{};
    for (int i = 0; i < Pop_Size; i++)
      for (int j = 0; j < Dim; j++)
        mean[j] += q * positions[i][j];
    std::stable_sort(index.begin(), index.end(),
                     [&](auto x, auto y) { return fit[x] > fit[y]; });
    if (fit[index.back()] < min_value) {
      min_value = fit[index.back()];
    }
    if constexpr (Dim > Pop_Size)
      std::generate_n(tmp_P.begin(), Pop_Size - 1, [&] { return sr1() * k; });
    if constexpr (Dim <= Pop_Size) {
      for (int i = 0; i < Pop_Size - 1; i++) {
        auto lo = index[i];
        for (int j = 0; j < Dim; j++) {
          win[j] = (sr2() % (Pop_Size - i - 1)) + i + 1;
          r1[j] = sr3() * k;
          r2[j] = sr4() * k;
          r3[j] = sr5() * k * eps;
        }
        for (int j = 0; j < Dim; j++) {
          auto w = index[win[j]];
          vec[lo][j] = r1[j] * vec[lo][j] +
                       r2[j] * (positions[w][j] - positions[lo][j]) +
                       r3[j] * (mean[j] - positions[lo][j]);
        }
        for (int j = 0; j < Dim; j++)
          positions[lo][j] = std::clamp(positions[lo][j] + vec[lo][j], l, r);
        fit[lo] = function(positions[lo].data());
      }
    } else {
      for (int i = 0; i < Pop_Size - 1; i++)
        if (tmp_P[i] < P[i]) {
          auto lo = index[i];
          for (int j = 0; j < Dim; j++) {
            win[j] = (sr2() % (Pop_Size - i - 1)) + i + 1;
            r1[j] = sr3() * k;
            r2[j] = sr4() * k;
            r3[j] = sr5() * k * eps;
          }
          for (int j = 0; j < Dim; j++) {
            auto w = index[win[j]];
            vec[lo][j] = r1[j] * vec[lo][j] +
                         r2[j] * (positions[w][j] - positions[lo][j]) +
                         r3[j] * (mean[j] - positions[lo][j]);
          }
          for (int j = 0; j < Dim; j++)
            positions[lo][j] = std::clamp(positions[lo][j] + vec[lo][j], l, r);
          fit[lo] = function(positions[lo].data());
        }
    }
    if constexpr (Dim <= Pop_Size) {
      fes += Pop_Size - 1;
    } else {
      for (int i = 0; i < Pop_Size - 1; i++)
        if (tmp_P[i] < P[i])
          fes++;
    }
    if constexpr (Memory_Flag) {
      slpso_convergence_curve[2 * _] = T(fes);
      slpso_convergence_curve[2 * _ + 1] =
          std::min(min_value, *std::min_element(fit.begin(), fit.end()));
    }
  }
  auto i=std::min_element(fit.begin(),fit.end())-fit.begin();
  if constexpr (Memory_Flag)
    return std::make_tuple(positions[i], fit[i],
                           std::move(slpso_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(positions[i], fit[i]);
}

template <int Dim, int Pop_Size, int Max, bool Memory_Flag = false, typename F,
          std::floating_point T, typename Parameter_Type = slpso_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           slpso_parameter_concept<Parameter_Type, T>
inline auto slpso(F &&function, T l, T r,
                  const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return slpso_optimize<Dim, Pop_Size, Max, Memory_Flag>(pop, function, l, r,
                                                         pt);
}

} // namespace sevobench
