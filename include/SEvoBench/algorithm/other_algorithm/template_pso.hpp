#pragma once
#include "../../common/common_concept.hpp"
#include "../../common/tool.hpp"

namespace sevobench::other_algorithm::pso_detail {

template <int Dim, std::floating_point T>
inline auto hypersphere_sample(T *first, T r, T *out) noexcept {
  auto sum = std::sqrt(std::inner_product(out, out + Dim, out, T(0)));
  for (int i = 0; i < Dim; i++)
    out[i] = first[i] + r * out[i] / sum;
  return out + Dim;
}

} // namespace sevobench::other_algorithm::pso_detail

namespace sevobench::other_algorithm {

template <std::floating_point T> struct spso2007_parameter {
  T w = T(1) / T(2 * std::numbers::ln2_v<T>);
  T c = T(0.5) + std::numbers::ln2_v<T>;
};
template <std::floating_point T> struct spso2011_parameter {
  T w = T(1) / T(2 * std::numbers::ln2_v<T>);
  T c = T(0.5) + std::numbers::ln2_v<T>;
};

template <typename P, typename T>
concept spso_parameter_concept = requires(const P &p) {
  requires std::same_as<T, std::remove_cvref_t<decltype(p.w)>>;
  requires std::same_as<T, std::remove_cvref_t<decltype(p.c)>>;
};

template <bool Version, int Dim, int Pop_Size, int Max, bool Memory_Flag,
          typename G, typename F, std::floating_point T,
          typename Parameter_Type = std::conditional_t<
              Version, spso2007_parameter<T>, spso2011_parameter<T>>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           spso_parameter_concept<Parameter_Type, T> &&
           algorithm_positions_concept<F, G, T>
inline auto
template_pso_optimize(G &&positions, F &&f, T left_bound, T right_bound,
                      const Parameter_Type &pt = Parameter_Type()) noexcept {
  const T w = pt.w;
  const T c = pt.c;
  constexpr T k = T(1) / T(0x7FFF);
  constexpr int K = Version ? 3 : 4;
  std::vector<std::array<T, Dim>> vec(Pop_Size);
  std::array<T, Pop_Size> pbest_fit;
  std::vector<std::array<int, K>> topology_index(Pop_Size);
  tool::curve_t<Memory_Flag, T, Max> template_pso_convergence_curve;
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<T> dis(left_bound, right_bound);
  auto pbest(positions);
  std::transform(positions.begin(), positions.begin() + Pop_Size,
                 pbest_fit.begin(), [&](auto &x) { return f(x.data()); });
  auto global_best_index =
      std::min_element(pbest_fit.begin(), pbest_fit.end()) - pbest_fit.begin();
  auto min_value = pbest_fit[global_best_index];
  tool::simple_rand sr1(rd());
  tool::simple_rand sr2(rd());
  tool::simple_rand sr3(rd());
  for (int i = 0; i < Pop_Size; i++) {
    if constexpr (Version) {
      for (int j = 0; j < K; j++)
        topology_index[i][j] = (i + j - 1) % Pop_Size;
      topology_index[0][0] = (Pop_Size - 1);
    } else {
      topology_index[i][0] = i;
      for (int j = 1; j < K; j++)
        topology_index[i][j] = sr1() % Pop_Size;
    }
  }
  for (int i = 0; i < Pop_Size; i++)
    for (int j = 0; j < Dim; j++) {
      if constexpr (Version) {
        vec[i][j] = (dis(gen) - positions[i][j]) * T(0.5);
      } else {
        vec[i][j] = (dis(gen) - positions[i][j]);
      }
    }
  for (int _ = (Memory_Flag ? 0 : Pop_Size); _ < Max;
       _ += (Memory_Flag ? 1 : Pop_Size)) {
    for (int i = 0; i < Pop_Size; i++) {
      auto index = *std::min_element(
          topology_index[i].begin(), topology_index[i].end(),
          [&](auto i1, auto i2) { return pbest_fit[i1] < pbest_fit[i2]; });
      if constexpr (Version) {
        for (int j = 0; j < Dim; j++)
          vec[i][j] =
              index == i
                  ? w * vec[i][j] +
                        c * sr2() * k * (pbest[i][j] - positions[i][j])
                  : w * vec[i][j] +
                        c * sr2() * k * (pbest[i][j] - positions[i][j]) +
                        c * sr3() * k * (pbest[index][j] - positions[i][j]);
      } else {
        std::array<T, Dim> gpos;
        std::array<T, Dim> xpos;
        T rsquare(0);
        for (int j = 0; j < Dim; j++) {
          T tmp =
              index == i
                  ? c * (pbest[i][j] - positions[i][j]) * T(0.5)
                  : c * (pbest[i][j] + pbest[index][j] - 2 * positions[i][j]) /
                        T(3);
          rsquare += tmp * tmp;
          gpos[j] = tmp + positions[i][j];
        }
        for (int j = 0; j < Dim; j++)
          xpos[j] = 2 * k * sr3() - 1;
        pso_detail::hypersphere_sample<Dim>(
            gpos.data(), std::sqrt(rsquare) * k * sr3(), xpos.data());
        for (int j = 0; j < Dim; j++)
          vec[i][j] = w * vec[i][j] + xpos[j] - positions[i][j];
      }
    }
    for (int i = 0; i < Pop_Size; i++) {
      for (int j = 0; j < Dim; j++) {
        positions[i][j] = positions[i][j] + vec[i][j];
        if (positions[i][j] < left_bound || positions[i][j] > right_bound) {
          positions[i][j] =
              std::clamp(positions[i][j], left_bound, right_bound);
          vec[i][j] = Version ? T(0) : -T(0.5) * vec[i][j];
        }
      }
      T tmp_fit = f(positions[i].data());
      if (tmp_fit < pbest_fit[i]) {
        pbest_fit[i] = tmp_fit;
        pbest[i] = positions[i];
      }
    }
    global_best_index = std::min_element(pbest_fit.begin(), pbest_fit.end()) -
                        pbest_fit.begin();
    if constexpr (!Version) {
      if (pbest_fit[global_best_index] >= min_value) {
        for (int i = 0; i < Pop_Size; i++) {
          topology_index[i][0] = i;
          for (int j = 1; j < K; j++)
            topology_index[i][j] = sr1() % Pop_Size;
        }
      } else {
        min_value = pbest_fit[global_best_index];
      }
    } else {
      if (pbest_fit[global_best_index] < min_value) {
        min_value = pbest_fit[global_best_index];
      }
    }
    if constexpr (Memory_Flag) {
      template_pso_convergence_curve[2 * _] = (_ + 1) * Pop_Size;
      template_pso_convergence_curve[2 * _ + 1] = min_value;
    }
  }
  if constexpr (Memory_Flag)
    return std::make_tuple(pbest[global_best_index], min_value,
                           std::move(template_pso_convergence_curve));
  if constexpr (!Memory_Flag)
    return std::make_pair(pbest[global_best_index], min_value);
}

template <int Dim, int Pop_Size = 30, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = spso2007_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           spso_parameter_concept<Parameter_Type, T>
inline auto spso2007(F &&function, T l, T r,
                     const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return template_pso_optimize<true, Dim, Pop_Size, Max, Memory_Flag>(
      pop, function, l, r, pt);
}
template <int Dim, int Pop_Size = 30, int Max = 1000 * Dim,
          bool Memory_Flag = false, typename F, std::floating_point T,
          typename Parameter_Type = spso2011_parameter<T>>
  requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
           spso_parameter_concept<Parameter_Type, T>
inline auto spso2011(F &&function, T l, T r,
                     const Parameter_Type &pt = Parameter_Type()) noexcept {
  auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
  return template_pso_optimize<false, Dim, Pop_Size, Max, Memory_Flag>(
      pop, function, l, r, pt);
}

} // namespace sevobench::other_algorithm
