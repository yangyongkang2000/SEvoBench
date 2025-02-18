#pragma once
#include "../../version2/vectormath_exp.h"
#include "../../version2/vectormath_trig.h"
#include "simd_type.hpp"
#include <numbers>
#include <numeric>

namespace sevobench::problem::simd {
template <int N, std::floating_point T>
inline void shift_func(const T *x, T *y, const T *o) noexcept {
  for (int i = 0; i < N; i++) {
    y[i] = x[i] - o[i];
  }
}

template <int N, std::floating_point T>
inline void rotate_func(T *x, const T *m) noexcept {
  std::array<T, N> tmp;
  std::copy_n(x, N, tmp.data());
  for (int i = 0; i < N; i++)
    x[i] = std::inner_product(tmp.data(), tmp.data() + N, m + i * N, T(0));
}

template <bool Using_SIMD, int N, int W, std::floating_point T>
inline auto sphere(const T *x) noexcept {
  return std::inner_product(x, x + N, x, T(0));
}

template <bool Using_SIMD, int N, int W, std::floating_point T>
inline auto rosenbrock(const T *x) noexcept {
  T sum(0);
  for (int i = 0; i < N - 1; i++)
    sum +=
        x[i] * x[i] + 100 * tool::Pow<2>(x[i + 1] + 1 - tool::Pow<2>(x[i] + 1));
  return sum;
}
template <bool Using_SIMD, int N, int W, std::floating_point T>
inline auto rastrigin(const T *x) noexcept {
  if constexpr (!Using_SIMD) {
    return 10 * N + std::accumulate(x, x + N, T(0), [](auto l, auto r) {
             return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
           });
  } else {
    vsimd_type<T, W> v(0);
    vsimd_type<T, W> w(0);
    T sum(10 * N);
    constexpr int RN = N & (-W);
    for (int i = 0; i < RN; i += W) {
      vsimd_type<T, W> tmp;
      tmp.load(x + i);
      v += tmp * tmp;
      w -= 10 * cospi(2 * tmp);
    }
    if constexpr (RN < N) {
      vsimd_type<T, W> tmp;
      tmp.load_partial(N - RN, x + RN);
      v += tmp * tmp;
      w -= 10 * cospi(2 * tmp);
    }
    sum += horizontal_add(v + w);
    if constexpr (RN < N) {
      sum += 10 * (W - (N - RN));
    }
    return sum;
  }
}

template <bool Using_SIMD, int N, int W, std::floating_point T>
inline auto ackley(const T *x) noexcept {
  if constexpr (!Using_SIMD) {
    T sum(T(0));
    T sqsum(T(0));
    for (int i = 0; i < N; i++) {
      sum += std::cos(2 * std::numbers::pi_v<T> * x[i]);
      sqsum += x[i] * x[i];
    }
    return 20 + std::numbers::e_v<T> - std::exp(sum / N) -
           20 * std::exp(T(-0.2) * std::sqrt(sqsum / N));
  } else {
    vsimd_type<T, W> sum(0);
    vsimd_type<T, W> sqsum(0);
    constexpr int RN = N & (-W);
    for (int i = 0; i < RN; i += W) {
      vsimd_type<T, W> tmp;
      tmp.load(x + i);
      sum += cospi(2 * tmp);
      sqsum += tmp * tmp;
    }
    if constexpr (RN < N) {
      vsimd_type<T, W> tmp;
      tmp.load_partial(N - RN, x + RN);
      sum += cospi(2 * tmp);
      sqsum += tmp * tmp;
    }
    T sq = horizontal_add(sqsum);
    T s = horizontal_add(sum);
    if constexpr (RN < N) {
      s -= (W - (N - RN));
    }
    return 20 + std::numbers::e_v<T> - std::exp(s / N) -
           20 * std::exp(T(-0.2) * std::sqrt(sq / T(N)));
  }
}

template <bool Using_SIMD, int N, int W, std::floating_point T>
inline auto elliptic(const T *x) noexcept {
  if constexpr (!Using_SIMD) {
    T sum(0);
    for (int i = 0; i < N; i++) {
      sum += std::pow(T(10), T(6 * i) / T(N - 1)) * x[i] * x[i];
    }
    return sum;
  } else {
    constexpr auto ite = [] {
      std::array<T, W> tmp{};
      std::iota(tmp.begin(), tmp.end(), T(0));
      return tmp;
    }();
    constexpr int RN = N & (-W);
    vsimd_type<T, W> v(0);
    for (int i = 0; i < RN; i += W) {
      vsimd_type<T, W> tmp;
      vsimd_type<T, W> tmp1;
      tmp1.load(ite.data());
      tmp1 += T(i);
      tmp.load(x + i);
      v += square(tmp) * pow(vsimd_type<T, W>(10), 6 * tmp1 / T(N - 1));
    }
    if constexpr (RN < N) {
      vsimd_type<T, W> remain;
      vsimd_type<T, W> tmp1;
      remain.load_partial(N - RN, x + RN);
      tmp1.load(ite.data());
      tmp1 += RN;
      v += square(remain) * pow(vsimd_type<T, W>(10), 6 * tmp1 / T(N - 1));
    }
    return horizontal_add(v);
  }
}

template <bool Using_SIMD, int N, int W, std::floating_point T>
inline auto schwefel_1(const T *x) noexcept {
  T sum(0);
  T s1(0);
  for (int i = 0; i < N; i++) {
    s1 += x[i];
    sum += s1 * s1;
  }
  return sum;
}

} // namespace sevobench::problem::simd