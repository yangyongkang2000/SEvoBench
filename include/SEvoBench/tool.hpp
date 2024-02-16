#pragma once

#include "common_concept.hpp"
#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <iterator>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace sevobench::tool {

class simple_rand {
  unsigned int g_seed;

public:
  simple_rand(unsigned int _) : g_seed(_) {}
  int operator()() noexcept {
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
  }
};

template <int N, std::floating_point T> inline constexpr T Pow(T x) noexcept {
  if constexpr (N <= 0)
    return 1;
  if constexpr (N % 2 == 0) {
    auto t = Pow<N / 2>(x);
    return t * t;
  } else {
    return x * Pow<(N - 1)>(x);
  }
}

template <int M, typename F>
inline auto binary_find_iteration(F &&f, int fes1, int fes2) noexcept {
  while (fes2 - fes1 > 1) {
    auto fes3 = (fes1 + fes2) / 2;
    auto ite = f(fes3);
    if (ite != M) {
      if (ite < M)
        fes1 = fes3;
      else
        fes2 = fes3;
    } else
      return fes3;
  }
  return fes1;
}

template <int N, std::floating_point T, typename V>
inline auto orthogonal_basis(V &b) noexcept {

  std::random_device rd;
  std::normal_distribution<T> norm(T(0), T(1));
  std::default_random_engine gen(rd()); // TODO: seed ?
  T sp(0);

  for (int i = 0; i < N; i++) {
    /* sample from Gaussian. */
    for (int j = 0; j < N; j++)
      b[i][j] = norm(gen);
    /* substract projection of previous vectors */
    for (int j = i - 1; j >= 0; --j) {
      sp = T(0);
      for (int k = 0; k < N; ++k)
        sp += b[i][k] * b[j][k]; // scalar product.
      for (int k = 0; k < N; k++)
        b[i][k] -= sp * b[j][k];
    }
    sp = T(0);
    for (int k = 0; k < N; ++k)
      sp += b[i][k] * b[i][k];
    for (int k = 0; k < N; ++k)
      b[i][k] /= std::sqrt(sp);
  }
}

template <typename InputIt>
inline constexpr auto find_median(InputIt begin, InputIt end) noexcept {
  using T = typename std::iterator_traits<InputIt>::value_type;
  auto count = end - begin;
  if (count % 2) {
    return begin[count / 2];
  } else {
    auto right = begin[count / 2];
    auto left = begin[count / 2 - 1];
    return (right + left) / T(2);
  }
}
template <typename InputIt>
inline auto mean_std(InputIt beg, InputIt end) noexcept {
  using T = typename std::iterator_traits<InputIt>::value_type;
  auto len = end - beg;
  auto mean = std::accumulate(beg, end, T(0)) / len;
  auto var = std::sqrt(std::inner_product(
                 beg, end, beg, T(0), std::plus<T>(),
                 [=](auto l, auto r) { return (l - mean) * (r - mean); })) /
             std::sqrt(len - 1);
  std::sort(beg, end);
  return std::make_tuple(
      mean, var, beg[0], tool::find_median(beg, beg + len / 2),
      tool::find_median(beg, end),
      tool::find_median(beg + len / 2 + (len % 2), end), beg[len - 1]);
}

template <class RandomIt>
inline void kunth_shuffle(RandomIt first, RandomIt last,
                          unsigned int seed) noexcept {
  typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
  simple_rand sr(seed);
  for (diff_t i = last - first - 1; i > 0; --i)
    std::swap(first[i], first[sr() % (i + 1)]);
}

template <std::floating_point T>
inline auto box_muller(T m, T st, T u, T v) noexcept {
  return st * std::sqrt(T(-2) * std::log(u)) *
             std::sin(2 * std::numbers::pi_v<T> * v) +
         m;
}
template <std::floating_point T>
inline auto cauchy_dis(T a, T b, T y) noexcept {
  return a + b * std::tan((y - T(0.5)) * std::numbers::pi_v<T>);
}

template <std::floating_point T, int Max>
class curve_vector : public std::vector<T> {
public:
  curve_vector() : std::vector<T>(2 * Max){};
};

template <bool b, std::floating_point T, int Max>
using curve_t = std::conditional_t<b, curve_vector<T, Max>, void *>;

template <std::floating_point T, int Dim, int Pop_Size>
inline auto random_generate_position(T l, T r) noexcept {
  std::vector<std::array<T, Dim>> positions(Pop_Size);
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<T> dis(l, r);
  std::for_each_n(positions.begin(), Pop_Size, [&](auto &x) {
    std::generate_n(x.begin(), Dim, [&] { return dis(gen); });
  });
  return positions;
}

template <std::unsigned_integral T>
constexpr T alg_hash(const char *s) noexcept {
  T hv(0);
  while (*s != '\0')
    hv = T(131 * hv) + T(*s++);
  return hv;
}

} // namespace sevobench::tool
