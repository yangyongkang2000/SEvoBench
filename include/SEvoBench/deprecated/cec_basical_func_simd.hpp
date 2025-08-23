#pragma once

#include "simd_type.hpp"
#include "tool.hpp"
#include "version2/vectormath_exp.h"
#include "version2/vectormath_trig.h"
namespace sevobench::ieee_cec_set::CEC_Base_Func_Simd {

template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline void shift_func(const T *x, T *y, const T *o) noexcept {
  if constexpr (simd_id() < 0) {
    for (int i = 0; i < N; i++) {
      y[i] = x[i] - o[i];
    }
  } else {
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> v1;
      simd_type_detail::simd_type<T> v2;
      v1.load(x + i);
      v2.load(o + i);
      (v1 - v2).store(y + i);
    }
  }
}

template <int N, typename T>
  requires simd_type_detail::simd_dim<T, N>
auto simd_inner_product(const T *x, const T *y) {
  if constexpr (simd_id() < 0) {
    return std::inner_product(x, x + N, y, T(0));
  } else {
    simd_type_detail::simd_type<T> v(0);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      simd_type_detail::simd_type<T> tmp1;
      tmp.load(x + i);
      tmp1.load(y + i);
      v += tmp * tmp1;
    }
  }
}

template <int N, std::floating_point T, typename M>
  requires simd_type_detail::simd_dim<T, N>
inline void rotate_func(T *x, const M &m) noexcept {
  std::array<T, N> tmp;
  std::copy_n(x, N, tmp.data());
  for (int i = 0; i < N; i++)
    x[i] = simd_inner_product<N>(tmp.data(), m[i].data(), T(0));
}

template <int N, std::floating_point T>
inline auto sphere(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    return std::inner_product(x, x + N, x, T(0));
  } else {
    simd_type_detail::simd_type<T> v(0);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      tmp.load(x + i);
      v += tmp * tmp;
    }
    return horizontal_add(v);
  }
}

template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto schwefel(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    return std::abs(*std::max_element(
        x, x + N, [](auto l, auto r) { return std::abs(l) < std::abs(r); }));
  } else {
    simd_type_detail::simd_type<T> v(0);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      tmp.load(x + i);
      v = max(v, abs(tmp));
    }
    return horizontal_max(v);
  }
}
template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto rosenbrock(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    T sum(0);
    for (int i = 0; i < N - 1; i++)
      sum += (x[i] - 1) * (x[i] - 1) +
             100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]);
    return sum;
  } else {
    simd_type_detail::simd_type<T> sum(0);
    for (int i = 0; i < N - simd_type_detail::simd_width<T>;
         i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      simd_type_detail::simd_type<T> tmp1;
      tmp.load(x + i);
      tmp1.load(x + i + 1);
      sum += square(tmp - 1) + 100 * square(tmp1 - tmp * tmp);
    }
    T s(0);
    for (int i = N - simd_type_detail::simd_width<T>; i < N - 1; i++)
      s += (x[i] - 1) * (x[i] - 1) +
           100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]);
    return horizontal_add(sum) + s;
  }
}
template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto rastrigin(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    return 10 * N + std::accumulate(x, x + N, T(0), [](auto l, auto r) {
             return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
           });
  } else {
    simd_type_detail::simd_type<T> v(0);
    simd_type_detail::simd_type<T> w(0);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      tmp.load(x + i);
      v += tmp * tmp;
      w -= 10 * cospi(2 * tmp);
    }
    return 10 * N + horizontal_add(v + w);
  }
}

template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto griewank(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    T s(0);
    T p(1);
    for (int i = 0; i < N; i++) {
      s += x[i] * x[i];
      p *= std::cos(x[i] / std::sqrt(T(1 + i)));
    }
    return 1 + s / T(4000) - p;
  } else {
    simd_type_detail::simd_type<T> s(0);
    simd_type_detail::simd_type<T> p(1);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      tmp.load(x + i);
      simd_type_detail::simd_type<T> tmp1;
      T data[simd_type_detail::simd_width<T>];
      s += tmp * tmp;
      std::iota(data, data + simd_type_detail::simd_width<T>, i + 1);
      tmp1.load(data);
      p *= cos(tmp / sqrt(tmp1));
    }
    T p1(1);
    for (int i = 0; i < simd_type_detail::simd_width<T>; i++)
      p1 *= p[i];
    return 1 - p1 + horizontal_add(s) / T(4000);
  }
}

template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto ackley(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    T sum(T(0));
    T sqsum(T(0));
    for (int i = 0; i < N; i++) {
      sum += std::cos(2 * std::numbers::pi_v<T> * x[i]);
      sqsum += x[i] * x[i];
    }
    return 20 + std::numbers::e_v<T> - std::exp(sum / N) -
           20 * std::exp(T(-0.2) * std::sqrt(sqsum / N));
  } else {
    simd_type_detail::simd_type<T> sum(0);
    simd_type_detail::simd_type<T> sqsum(0);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      tmp.load(x + i);
      sum += cospi(2 * tmp);
      sqsum += tmp * tmp;
    }
    return 20 + std::numbers::e_v<T> - std::exp(horizontal_add(sum) / N) -
           20 * std::exp(T(-0.2) * std::sqrt(horizontal_add(sqsum) / T(N)));
  }
}

template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto elliptic(const T *x) noexcept {
  if constexpr (simd_id() < 0) {
    T sum(0);
    for (int i = 0; i < N; i++) {
      sum += std::pow(T(10), T(6 * i) / T(N - 1)) * x[i] * x[i];
    }
    return sum;
  } else {
    simd_type_detail::simd_type<T> sum(0);
    T y[N];
    y[0] = 1;
    auto ratio = std::pow(1000000, T(1) / T(N - 1));
    for (int i = 1; i < N; i++)
      y[i] = y[i - 1] * ratio;
    simd_type_detail::simd_type<T> v(0);
    for (int i = 0; i < N; i += simd_type_detail::simd_width<T>) {
      simd_type_detail::simd_type<T> tmp;
      simd_type_detail::simd_type<T> tmp1;
      tmp.load(x + i);
      tmp1.load(y + i);
      v += tmp * tmp1;
    }
    return horizontal_add(v);
  }
}

template <int N, std::floating_point T>
  requires simd_type_detail::simd_dim<T, N>
inline auto schwefel_1(const T *x) noexcept {
  std::array<T, N> tmp;
  std::partial_sum(x, x + N, tmp.data());
  return sphere<N>(tmp.data());
}

} // namespace sevobench::ieee_cec_set::CEC_Base_Func_Simd
