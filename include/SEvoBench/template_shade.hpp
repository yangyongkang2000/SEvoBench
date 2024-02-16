#pragma once

#include "tool.hpp"

namespace sevobench::shade_detail {
template <int N_Max, int N_Min>
inline auto lshade_ite(const int fes_max) noexcept {
  int fes(N_Max);
  int n(N_Max);
  int i(0);
  for (; fes < fes_max; ++i) {
    fes += n;
    n = std::max(static_cast<int>(std::round(
                     N_Max + (N_Min - N_Max) / float(fes_max) * float(fes))),
                 N_Min);
  }
  return i;
}

template <int N_Max, int N_Min, int Ite>
static const int iteration_maxfes = tool::binary_find_iteration<Ite>(
    lshade_ite<N_Max, N_Min>, N_Min *Ite, N_Max *Ite);

} // namespace sevobench::shade_detail

namespace sevobench {

template <std::floating_point T> struct shade_parameter {
  T r_Arc{2};
  T p = T(0.1);
  int H{100};
};
template <std::floating_point T> struct lshade_parameter {
  T r_Arc{2.6};
  T p = T(0.11);
  int H{6};
};

template <typename P, typename T>
concept template_shade_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.r_Arc)>>;
  requires std::same_as<T, std::remove_cvref_t<decltype(p.p)>>;
  requires std::same_as<int, std::remove_cvref_t<decltype(p.H)>>;
};

template <bool b, int Dim, int Pop_Size, int Max, bool Memory_Flag, typename G,
          typename F, std::floating_point T,
          typename Parameter_Type =
              std::conditional_t<b, shade_parameter<T>, lshade_parameter<T>>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           template_shade_parameter_concept<Parameter_Type, T> &&
           algorithm_positions_concept<F, G, T> && (Pop_Size >= 3)
inline auto template_shade_optimize(
    G &&positions, F &&function, T l, T r,
    const Parameter_Type &pt = Parameter_Type()) {
  const int LP = pt.H;
  constexpr int N_Min = 4;
  int M_Gen = Pop_Size;
  int A_Size = static_cast<int>(pt.r_Arc * Pop_Size);
  constexpr T k = T(1) / T(0x7fff);
  std::vector<T> MCR(LP, T(0.5));
  std::vector<T> MF(LP, T(0.5));
  std::vector<T> f_archive(Pop_Size);
  std::vector<T> cr_archive(Pop_Size);
  std::vector<T> delta_f(Pop_Size);
  std::vector<std::array<T, Dim>> Archive(A_Size);
  int Archive_Size = 0;
  int replace_index = 0;
  int index_counter = 0;
  std::vector<std::array<T, Dim>> tmp(M_Gen);
  tool::curve_t<Memory_Flag, T, Max> shade_convergence_curve;
  std::vector<T> fit(Pop_Size);
  std::vector<T> tmp_fit(Pop_Size);
  std::array<T, Dim> random_c;
  std::vector<int> index(Pop_Size);
  std::random_device rd;
  tool::simple_rand sr1(rd());
  tool::simple_rand sr2(rd());
  tool::simple_rand sr3(rd());
  tool::simple_rand sr4(rd());
  tool::simple_rand sr5(rd());
  std::iota(index.begin(), index.end(), 0);
  std::transform(positions.begin(), positions.begin() + Pop_Size, fit.begin(),
                 [&](auto &x) { return function(x.data()); });
  auto gen_f = [&](auto x, auto &g) {
    if (!std::isfinite(x))
      return T(1);
    auto z = tool::cauchy_dis(x, T(0.1), g() * k);
    constexpr int sample_max = 20;
    for (int i = 0; z <= T(0); i++) {
      z = tool::cauchy_dis(x, T(0.1), g() * k);
      if (i > sample_max)
        return T(1);
    }
    return z < 1 ? z : T(1);
  };
  int fes(Pop_Size);
  for (int _ = 0; (Memory_Flag ? _ : fes) < Max; _++) {
    int P =
        static_cast<int>(b ? Pop_Size * pt.p : std::max(M_Gen * pt.p, T(2)));
    auto f_first = f_archive.begin();
    auto cr_first = cr_archive.begin();
    auto df_first = delta_f.begin();
    std::nth_element(index.begin(), index.begin() + P, index.begin() + M_Gen,
                     [&](auto x, auto y) { return fit[x] < fit[y]; });
    for (int i = 0; i < M_Gen; i++) {
      auto rp = sr1() % LP;
      auto cR =
          (MCR[rp] == T(-1) ? T(0)
                            : std::clamp(tool::box_muller(MCR[rp], T(0.1),
                                                          sr1() * k, sr2() * k),
                                         T(0), T(1)));
      auto cF = gen_f(MF[rp], sr3);
      int r0, r1, r2;
      r0 = index[P == 0 ? 0 : sr4() % P];
      do {
        r1 = sr4() % M_Gen;
      } while (r1 == i);
      do {
        r2 = sr4() % (M_Gen + Archive_Size);
      } while (r2 == i || r2 == r1);
      auto j_rand = sr4() % Dim;
      std::generate_n(random_c.begin(), Dim, [&] { return sr5() * k; });
      auto &vr2 =
          r2 >= M_Gen ? Archive[r2 - M_Gen] : positions[b ? r2 : index[r2]];
      r1 = b ? r1 : index[r1];
      auto i1 = b ? i : index[i];
      for (int j = 0; j < Dim; j++) {
        auto x = (j == j_rand || random_c[j] < cR)
                     ? positions[i1][j] +
                           cF * (positions[r0][j] - positions[i1][j]) +
                           cF * (positions[r1][j] - vr2[j])
                     : positions[i1][j];
        if (x < l)
          x = (positions[i1][j] + l) * T(0.5);
        if (x > r)
          x = (positions[i1][j] + r) * T(0.5);
        tmp[i][j] = x;
      }
      tmp_fit[i] = function(tmp[i].data());
      if (tmp_fit[i] < fit[i1]) {
        *f_first++ = cF;
        *cr_first++ = cR;
        *df_first++ = fit[i1] - tmp_fit[i];
        Archive[replace_index++] = positions[i1];
        replace_index %= A_Size;
        Archive_Size = Archive_Size == A_Size ? A_Size : Archive_Size + 1;
      }
    }
    for (int i = 0; i < M_Gen; i++)
      if (tmp_fit[i] < fit[b ? i : index[i]]) {
        fit[b ? i : index[i]] = tmp_fit[i];
        positions[b ? i : index[i]] = tmp[i];
      }
    if (f_first != f_archive.begin()) {
      auto scr2 = std::inner_product(cr_archive.begin(), cr_first,
                                     delta_f.begin(), T(0), std::plus<T>(),
                                     [](auto x, auto y) { return x * x * y; });
      auto scr1 = std::inner_product(cr_archive.begin(), cr_first,
                                     delta_f.begin(), T(0));
      MCR[index_counter] = (scr1 == T(0) ? T(-1) : scr2 / scr1);
      MF[index_counter] =
          std::inner_product(f_archive.begin(), f_first, delta_f.begin(), T(0),
                             std::plus<T>(),
                             [](auto x, auto y) { return x * x * y; }) /
          std::inner_product(f_archive.begin(), f_first, delta_f.begin(), T(0));
      index_counter = (index_counter + 1) % LP;
    }
    fes += M_Gen;
    if constexpr (!b) {
      int tmp_M_Gen = std::max(
          static_cast<int>(std::round(
              T(N_Min - Pop_Size) * T(fes) /
                  T(Memory_Flag
                        ? shade_detail::iteration_maxfes<Pop_Size, N_Min, Max>
                        : Max) +
              T(Pop_Size))),
          N_Min);
      if (tmp_M_Gen < M_Gen) {
        std::nth_element(index.begin(), index.begin() + tmp_M_Gen,
                         index.begin() + M_Gen,
                         [&](auto x, auto y) { return fit[x] < fit[y]; });
        A_Size = static_cast<int>(pt.r_Arc * tmp_M_Gen);
        std::reverse(Archive.begin(), Archive.begin() + Archive_Size);
        Archive_Size = std::min(A_Size, Archive_Size);
        replace_index %= Archive_Size;
        M_Gen = tmp_M_Gen;
      }
    }
    if constexpr (Memory_Flag) {
      shade_convergence_curve[2 * _] = fes;
      shade_convergence_curve[2 * _ + 1] =
          *std::min_element(fit.begin(), fit.end());
    }
  }
  auto i = std::min_element(fit.begin(), fit.end()) - fit.begin();
  if constexpr (Memory_Flag)
    return std::make_tuple(positions[i], fit[i],
                           std::move(shade_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(positions[i], fit[i]);
}

template <int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = shade_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           template_shade_parameter_concept<Parameter_Type, T> &&
           (Pop_Size >= 3)
inline auto shade(F &&function, T l, T r,
                  const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return template_shade_optimize<true, Dim, Pop_Size, Max, Memory_Flag>(
      pop, function, l, r, pt);
}

template <int Dim, int Pop_Size = 18 * Dim, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = lshade_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           template_shade_parameter_concept<Parameter_Type, T> &&
           (Pop_Size >= 3)
inline auto lshade(F &&function, T l, T r,
                   const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return template_shade_optimize<false, Dim, Pop_Size, Max, Memory_Flag>(
      pop, function, l, r, pt);
}

} // namespace sevobench
