#pragma once
#include "../../common/tool.hpp"

namespace sevobench::problem::cec_detail {

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto zakharov_func(std::span<T, M> y) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? y.size() : M);
  T sum1(0);
  T sum2(0);
  for (int i = 0; i < N; i++) {
    sum1 += y[i] * y[i];
    sum2 += T(0.5) * (i + 1) * y[i];
  }
  return sum1 + tool::Pow<2>(sum2) + tool::Pow<4>(sum2);
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto sum_diff_pow_func(std::span<T, M> y) noexcept {
  T sum = 0;
  const auto N = static_cast<int>((M == std::dynamic_extent) ? y.size() : M);
  for (int i = 0; i < N; i++) {
    T newv = std::pow(std::abs(y[i]), (i + 1));
    sum += newv;
  }
  return sum;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto rosenbrock_func(std::span<T, M> t) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? t.size() : M);
  T sum(0);
  for (int i = 0; i < N; i++)
    t[i] = T(1) + T(0.02048) * t[i];
  for (int i = 0; i < N - 1; i++)
    sum +=
        tool::Pow<2>(t[i] - 1) + 100 * tool::Pow<2>((t[i + 1] - t[i] * t[i]));
  return sum;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto escaffer6_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  T sum(0);
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

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto rastrigin_func(std::span<T, M> _x) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? _x.size() : M);
  return 10 * N +
         std::accumulate(_x.begin(), _x.end(), T(0), [](auto l, auto r) {
           r *= T(0.0512);
           return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
         });
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto levy_func(std::span<T, M> x) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? x.size() : M);
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

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto bent_cigar_func(std::span<T, M> x) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? x.size() : M);
  T sum(0);
  for (int i = 1; i < N; i++)
    sum += x[i] * x[i];
  return x[0] * x[0] + 1000000 * sum;
}
template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto hgbat_func(std::span<T, M> x) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? x.size() : M);
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

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto katsuura_func(std::span<T, M> z) noexcept {
  constexpr T b[] = {T(2),         T(4),          T(8),          T(16),
                     T(32),        T(64),         T(128),        T(256),
                     T(512),       T(1024),       T(2048),       T(4096),
                     T(8192),      T(16384),      T(32768),      T(65536),
                     T(131072),    T(262144),     T(524288),     T(1048576),
                     T(2097152),   T(4194304),    T(8388608),    T(16777216),
                     T(33554432),  T(67108864),   T(134217728),  T(268435456),
                     T(536870912), T(1073741824), T(2147483648), T(4294967296)};
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  T sum(1);
  T tmp3 = std::pow(T(N), T(1.2));
  for (int i = 0; i < N; i++)
    z[i] *= T(0.05);
  for (int i = 0; i < N; i++) {
    T temp(0);
    for (int j = 1; j <= 32; j++) {
      T tmp1 = b[j - 1];
      T tmp2 = tmp1 * z[i];
      temp += std::abs(tmp2 - std::floor(tmp2 + T(0.5))) / tmp1;
    }
    sum *= std::pow(T(1) + (i + 1) * temp, T(10) / tmp3);
  }
  T tmp1 = T(10) / T(N * N);
  sum = sum * tmp1 - tmp1;
  return sum;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto ackley_func(std::span<T, M> _x) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? _x.size() : M);
  T sum(T(0));
  T sqsum(T(0));
  for (int i = 0; i < N; i++) {
    sum += std::cos(2 * std::numbers::pi_v<T> * _x[i]);
    sqsum += _x[i] * _x[i];
  }
  return 20 + std::numbers::e_v<T> - std::exp(sum / N) -
         20 * std::exp(T(-0.2) * std::sqrt(sqsum / N));
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto schwefel_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  T sum(0);
  for (int i = 0; i < N; i++)
    z[i] *= T(10);
  for (int i = 0; i < N; i++) {
    T x(z[i] + T(420.9687462275036));
    if (x > T(500)) {
      sum -= (T(500) - std::fmod(x, T(500))) *
             std::sin(std::sqrt(std::abs(500 - std::fmod(x, T(500)))));
      T tmp = (x - T(500)) / T(100);
      sum += tmp * tmp / N;
    } else if (x < T(-500)) {
      sum -= (T(-500) + std::fmod(std::abs(x), T(500))) *
             std::sin(
                 std::sqrt(std::abs(T(500) - std::fmod(std::abs(x), T(500)))));
      T tmp = (x + T(500)) / T(100);
      sum += tmp * tmp / N;
    } else
      sum -= x * std::sin(std::sqrt(std::abs(x)));
  }
  sum += T(418.9828872724338) * N;
  return sum;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto schwefel_F7_func(std::span<T, M> y) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? y.size() : M);
  T sum(0);
  for (int i = 0; i < N - 1; i++) {
    auto x = std::sqrt(y[i] * y[i] + y[i + 1] * y[i + 1]);
    auto tmp = std::sin(T(50) * std::pow(x, T(0.2)));
    sum += std::sqrt(x) * (1 + tmp * tmp);
  }
  sum = sum * sum / T((N - 1) * (N - 1));
  return sum;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto happycat_func(std::span<T, M> y) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? y.size() : M);
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

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto grie_rosen_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
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
template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto griewank_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
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
template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto ellips_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  T sum(0);
  for (int i = 0; i < N; i++) {
    sum += std::pow(T(10), T(6 * i) / T(N - 1)) * z[i] * z[i];
  }
  return sum;
}
template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto discus_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  T sum(1000000 * z[0] * z[0]);
  for (int i = 1; i < N; i++) {
    sum += z[i] * z[i];
  }
  return sum;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto step_rastrigin_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  for (int i = 0; i < N; i++)
    z[i] = T(0.0512) * z[i];
  return 10 * N +
         std::accumulate(z.data(), z.data() + N, T(0), [](auto l, auto r) {
           r = std::abs(r) >= T(0.5) ? T(0.5) * std::round(2 * r) : r;
           return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
         });
}

template <int N, std::floating_point T>
inline auto bi_rastrigin_func(std::span<const T> z, std::span<const T> o,
                              std::span<const T> m) noexcept {
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
    y[i] =
        std::inner_product(tmp.data(), tmp.data() + N, m.data() + i * N, T(0));
  T sum3(N);
  for (int i = 0; i < N; i++)
    sum3 -= std::cos(2 * std::numbers::pi_v<T> * y[i]);
  return std::min(sum1, sum2) + 10 * sum3;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto bi_rastrigin_func(std::span<T, M> tmp) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? tmp.size() : M);
  constexpr T mu0(2.5);
  constexpr T d(1);
  const T s = T(1) - T(1) / (T(2) * std::sqrt(T(N + 20)) - T(8.2));
  const T mu1 = -std::sqrt((mu0 * mu0 - d) / s);
  for (int i = 0; i < N; i++)
    tmp[i] += mu0;
  T sum1(0);
  T sum2(N);
  for (int i = 0; i < N; i++) {
    sum1 += (tmp[i] - mu0) * (tmp[i] - mu0);
    sum2 += s * (tmp[i] - mu1) * (tmp[i] - mu1);
    tmp[i] -= mu0;
  }
  T sum3(N);
  for (int i = 0; i < N; i++)
    sum3 -= std::cos(2 * std::numbers::pi_v<T> * tmp[i]);
  return std::min(sum1, sum2) + 10 * sum3;
}

template <std::floating_point T, std::size_t M = std::dynamic_extent>
inline auto weierstrass_func(std::span<T, M> z) noexcept {
  const auto N = static_cast<int>((M == std::dynamic_extent) ? z.size() : M);
  constexpr T a[] = {T(1) / T(1),      T(1) / T(2),      T(1) / T(4),
                     T(1) / T(8),      T(1) / T(16),     T(1) / T(32),
                     T(1) / T(64),     T(1) / T(128),    T(1) / T(256),
                     T(1) / T(512),    T(1) / T(1024),   T(1) / T(2048),
                     T(1) / T(4096),   T(1) / T(8192),   T(1) / T(16384),
                     T(1) / T(32768),  T(1) / T(65536),  T(1) / T(131072),
                     T(1) / T(262144), T(1) / T(524288), T(1) / T(1048576)};
  constexpr T b[] = {
      tool::Pow<0, T>(3),  tool::Pow<1, T>(3),  tool::Pow<2, T>(3),
      tool::Pow<3, T>(3),  tool::Pow<4, T>(3),  tool::Pow<5, T>(3),
      tool::Pow<6, T>(3),  tool::Pow<7, T>(3),  tool::Pow<8, T>(3),
      tool::Pow<9, T>(3),  tool::Pow<10, T>(3), tool::Pow<11, T>(3),
      tool::Pow<12, T>(3), tool::Pow<13, T>(3), tool::Pow<14, T>(3),
      tool::Pow<15, T>(3), tool::Pow<16, T>(3), tool::Pow<17, T>(3),
      tool::Pow<18, T>(3), tool::Pow<19, T>(3), tool::Pow<20, T>(3)};
  constexpr int k_max = 20;
  T sum1 = 0;
  constexpr T sum2 = T(-2097151) / T(1048576);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j <= k_max; j++) {
      T x = T(0.5 / 100) * z[i];
      sum1 += a[j] * std::cos(2 * std::numbers::pi_v<T> * b[j] * (x + T(0.5)));
    }
  }
  sum1 -= N * sum2;
  return sum1;
}
} // namespace sevobench::problem::cec_detail
