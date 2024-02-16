#pragma once
#include "tool.hpp"

namespace sevobench::ieee_cec_set {
namespace CEC_Base_Func {

template <int N, std::floating_point T, typename M>
inline void sr_func(const T *x, T *y, const T *o, const M &m) noexcept {
  std::array<T, N> tmp;
  for (int i = 0; i < N; i++)
    tmp[i] = x[i] - o[i];
  for (int i = 0; i < N; i++)
    y[i] = std::inner_product(tmp.data(), tmp.data() + N, m[i].begin(), T(0));
}

template <int N, int Num, std::floating_point T>
inline auto cf_cal(const T *x, const T *o, T *fit, T *delta) noexcept {
  std::array<T, Num> w{};
  for (int i = 0; i < Num; i++) {
    for (int j = 0; j < N; j++)
      w[i] += tool::Pow<2>(x[j] - o[i * N + j]);
  }
  auto it = std::find(w.begin(), w.end(), T(0));
  if (it != w.end()) {
    return fit[it - w.begin()];
  }
  for (int i = 0; i < Num; i++)
    w[i] = std::exp(-w[i] / (2 * N * tool::Pow<2>(delta[i]))) / std::sqrt(w[i]);
  auto w_sum = std::accumulate(w.begin(), w.end(), T(0));
  return std::inner_product(w.begin(), w.end(), fit, T(0)) / w_sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto zakharov_func(T *y) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum1(0);
  T sum2(0);
  for (int i = 0; i < N; i++) {
    sum1 += y[i] * y[i];
    sum2 += T(0.5) * (i + 1) * y[i];
  }
  return sum1 + tool::Pow<2>(sum2) + tool::Pow<4>(sum2);
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto rosenbrock_func(T *t) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(0);
  for (int i = 0; i < N; i++)
    t[i] = T(1) + T(0.02048) * t[i];
  for (int i = 0; i < N - 1; i++)
    sum += (t[i] - 1) * (t[i] - 1) +
           100 * (t[i + 1] - t[i] * t[i]) * (t[i + 1] - t[i] * t[i]);
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto escaffer6_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(0);
  for (int i = 0; i < N; i++)
    z[i] = T(0.005) * z[i];
  for (int i = 0; i < N - 1; i++) {
    T temp1 = std::sin(std::sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
    temp1 = temp1 * temp1;
    T temp2 = T(1) + T(0.001) * (z[i] * z[i] + z[i + 1] * z[i + 1]);
    sum += T(0.5) + (temp1 - T(0.5)) / (temp2 * temp2);
  }
  T temp1 = std::sin(std::sqrt(z[N - 1] * z[N - 1] + z[0] * z[0]));
  temp1 = temp1 * temp1;
  T temp2 = T(1) + T(0.001) * (z[N - 1] * z[N - 1] + z[0] * z[0]);
  sum += T(0.5) + (temp1 - T(0.5)) / (temp2 * temp2);
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto rastrigin_func(T *_x) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  return 10 * N + std::accumulate(_x, _x + N, T(0), [](auto l, auto r) {
           r *= T(0.0512);
           return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
         });
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto levy_func(T *x) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  auto sum =
      tool::Pow<2>(std::sin(std::numbers::pi_v<T> * (1 + T(0.25) * (x[0]))));
  for (int i = 1; i < N - 1; i++) {
    auto w = 1 + T(0.25) * (x[i]);
    sum += (w - 1) * (w - 1) *
           (1 + 10 * tool::Pow<2>(std::sin(std::numbers::pi_v<T> * w + 1)));
  }
  auto w = 1 + T(0.25) * (x[N - 1]);
  sum += (w - 1) * (w - 1) *
         (1 + tool::Pow<2>(std::sin(2 * std::numbers::pi_v<T> * w)));
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto bent_cigar_func(T *x) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(x[0] * x[0]);
  for (int i = 1; i < N; i++)
    sum += 1000000 * x[i] * x[i];
  return sum;
}
template <int N, std::floating_point T>
  requires(N >= 0)
inline auto hgbat_func(T *x) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum1(0);
  T sum2(0);
  for (int i = 0; i < N; i++)
    x[i] *= T(0.05);
  for (int i = 0; i < N; i++) {
    auto z = x[i] - 1;
    sum1 += z * z;
    sum2 += z;
  }
  return std::sqrt(std::abs(sum1 * sum1 - sum2 * sum2)) +
         (T(0.5) * sum1 + sum2) / N + T(0.5);
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto katsuura_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(1);
  T tmp3 = std::pow(T(N), T(1.2));
  for (int i = 0; i < N; i++)
    z[i] *= T(0.05);
  for (int i = 0; i < N; i++) {
    T temp(0);
    for (int j = 1; j <= 32; j++) {
      T tmp1 = std::pow(T(2), T(j));
      T tmp2 = tmp1 * z[i];
      temp += std::abs(tmp2 - std::floor(tmp2 + T(0.5))) / tmp1;
    }
    sum *= std::pow(T(1) + (i + 1) * temp, T(10) / tmp3);
  }
  T tmp1 = T(10) / T(N * N);
  sum = sum * tmp1 - tmp1;
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto ackley_func(T *_x) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(T(0));
  T sqsum(T(0));
  for (int i = 0; i < N; i++) {
    sum += std::cos(2 * std::numbers::pi_v<T> * _x[i]);
    sqsum += _x[i] * _x[i];
  }
  return 20 + std::numbers::e_v<T> - std::exp(sum / N) -
         20 * std::exp(T(-0.2) * std::sqrt(sqsum / N));
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto schwefel_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(0);
  for (int i = 0; i < N; i++)
    z[i] *= T(10);
  for (int i = 0; i < N; i++) {
    T x(z[i] + T(420.9687462275036));
    if (x > 500) {
      sum -= (500 - std::fmod(x, T(500))) *
             std::sin(std::sqrt(std::abs(500 - std::fmod(x, T(500)))));
      T tmp = (x - 500) / T(100);
      sum += tmp * tmp / N;
    } else if (x < -500) {
      sum -= (-500 + std::fmod(std::abs(x), 500)) *
             std::sin(std::sqrt(std::abs(500 - std::fmod(std::abs(x), 500))));
      T tmp = (x + 500) / 100;
      sum += tmp * tmp / N;
    } else
      sum -= x * std::sin(std::sqrt(std::abs(x)));
  }
  sum += T(418.9828872724338) * N;
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto schwefel_F7_func(T *y) noexcept {
  if constexpr (N == 0 || N == 1) {
    return T(0);
  }
  T sum(0);
  for (int i = 0; i < N - 1; i++) {
    auto x = std::sqrt(y[i] * y[i] + y[i + 1] * y[i + 1]);
    auto tmp = std::sin(T(50) * std::pow(x, T(0.2)));
    sum += std::sqrt(x) * (1 + tmp);
  }
  sum = sum * sum / T((N - 1) * (N - 1));
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto happycat_func(T *y) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum1(0);
  T sum2(0);
  for (int i = 0; i < N; i++)
    y[i] *= T(0.05);
  for (int i = 0; i < N; i++) {
    auto z = y[i] - 1;
    sum1 += z * z;
    sum2 += z;
  }
  return std::pow(std::abs(sum1 - N), T(0.25)) + (T(0.5) * sum1 + sum2) / N +
         T(0.5);
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto grie_rosen_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(0);
  for (int i = 0; i < N; i++)
    z[i] *= T(0.05);
  z[0] += 1;
  for (int i = 0; i < N - 1; i++) {
    z[i + 1] += 1;
    auto tmp1 = z[i] * z[i] - z[i + 1];
    auto tmp2 = z[i] - 1;
    auto temp = T(100) * tmp1 * tmp1 + tmp2 * tmp2;
    sum += (temp * temp) / T(4000) - std::cos(temp) + T(1);
  }
  auto tmp1 = z[N - 1] * z[N - 1] - z[0];
  auto tmp2 = z[N - 1] - 1;
  auto temp = T(100) * tmp1 * tmp1 + tmp2 * tmp2;
  sum += (temp * temp) / T(4000) - std::cos(temp) + T(1);
  return sum;
}
template <int N, std::floating_point T>
  requires(N >= 0)
inline auto griewank_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T s(0);
  T p(1);
  for (int i = 0; i < N; i++)
    z[i] *= T(6);
  for (int i = 0; i < N; i++) {
    s += z[i] * z[i];
    p *= std::cos(z[i] / std::sqrt(T(1 + i)));
  }
  return 1 + s / T(4000) - p;
}
template <int N, std::floating_point T>
  requires(N >= 0)
inline auto ellips_func(T *z) noexcept {
  if constexpr (N == 0 || N == 1) {
    return N == 0 ? 0 : 1000000 * z[0];
  }
  T sum(0);
  for (int i = 0; i < N; i++) {
    sum += std::pow(T(10), T(6 * i) / T(N - 1)) * z[i] * z[i];
  }
  return sum;
}
template <int N, std::floating_point T>
  requires(N >= 0)
inline auto discus_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  T sum(1000000 * z[0] * z[0]);
  for (int i = 1; i < N; i++) {
    sum += z[i] * z[i];
  }
  return sum;
}

template <int N, std::floating_point T>
  requires(N >= 0)
inline auto step_rastrigin_func(T *z) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  for (int i = 0; i < N; i++)
    z[i] = T(0.0512) * z[i];
  return 10 * N + std::accumulate(z, z + N, T(0), [](auto l, auto r) {
           r = std::abs(r) >= T(0.5) ? T(0.5) * std::round(2 * r) : r;
           return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
         });
}

template <int N, std::floating_point T, typename M>
inline auto bi_rastrigin_func(const T *z, const T *o, const M &m) noexcept {
  if constexpr (N == 0) {
    return T(0);
  }
  constexpr T mu0(2.5);
  constexpr T d(1);
  T s = T(1) - T(1) / (T(2) * std::sqrt(T(N + 20)) - T(8.2));
  T mu1 = -std::sqrt((mu0 * mu0 - d) / s);
  std::array<T, N> tmp;
  std::array<T, N> y;
  for (int i = 0; i < N; i++)
    tmp[i] = T(0.2) * T(o[i] > 0 ? 1 : -1) * (z[i] - o[i]) + mu0;
  T sum1(0);
  T sum2(N);
  for (int i = 0; i < N; i++) {
    sum1 += (tmp[i] - mu0) * (tmp[i] - mu0);
    sum2 += s * (tmp[i] - mu1) * (tmp[i] - mu1);
    tmp[i] -= mu0;
  }
  for (int i = 0; i < N; i++)
    y[i] = std::inner_product(tmp.data(), tmp.data() + N, m[i].begin(), T(0));
  T sum3(N);
  for (int i = 0; i < N; i++)
    sum3 -= std::cos(2 * std::numbers::pi_v<T> * y[i]);
  return std::min(sum1, sum2) + 10 * sum3;
}

} // namespace CEC_Base_Func
} // namespace sevobench::ieee_cec_set
