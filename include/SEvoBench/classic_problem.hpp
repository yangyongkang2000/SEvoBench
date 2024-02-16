#pragma once

#include "tool.hpp"
namespace sevobench::yao_func {

namespace detail {

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
bool is_init_shift = false;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
bool is_init_rotate = false;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<T, Dim>> shift_vector_data(13);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> rotate_matrix_data(13);

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
inline void symmetry_breaking(T *x) noexcept {
  for (int i = 0; i < Dim; i++) {
    auto tmp = x[i] == T(0) ? 0 : std::log(std::abs(x[i]));
    x[i] = x[i] *
           std::exp(T(0.049) * (std::sin((x[i] > 0 ? T(10) : T(5.5)) * tmp) +
                                std::sin((x[i] > 0 ? T(7.9) : T(3.1)) * tmp)));
    if constexpr (Dim == 1) {
      x[i] = x[i] > 0 ? x[i] * pow(x[i], T(0.2) * std::sqrt(x[i])) : x[i];
      x[i] = std::sqrt(T(10)) * x[i];
    } else {
      x[i] = x[i] > 0
                 ? x[i] * pow(x[i], T(0.2 * i) / (Dim - 1) * std::sqrt(x[i]))
                 : x[i];
      x[i] = std::pow(10, T(0.5) * i / (Dim - 1)) * x[i];
    }
  }
}

} // namespace detail

template <bool b, int Dim, std::floating_point T>
  requires(Dim >= 1)
void init_shift_rotate() {
  if (!detail::is_init_shift<T, Dim>) {
    constexpr T num[] = {100, 10,      100, 100, 30, 100, T(1.28),
                         100, T(5.12), 32,  600, 50, 50};
    std::random_device rd;
    std::default_random_engine gen(rd());
    for (int i = 0; i < 13; i++) {
      std::uniform_real_distribution<T> dis(-num[i], num[i]);
      std::generate_n(detail::shift_vector_data<T, Dim>[i].begin(), Dim,
                      [&] { return dis(gen); });
    }
    detail::is_init_shift<T, Dim> = true;
  }
  if constexpr (b) {
    if (!detail::is_init_rotate<T, Dim>) {
      for (auto &_ : detail::rotate_matrix_data<T, Dim>)
        tool::orthogonal_basis<Dim, T>(_);
      detail::is_init_rotate<T, Dim> = true;
    }
  }
}
template <int Prob_Index, int Dim, std::floating_point T, bool b = false,
          bool asymmetry = false>
struct classic_problem;

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<1, Dim, T, b, asymmetry> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++)
      tmp[i] = x[i] - shift[i];
    auto _x = std::begin(tmp);
    if constexpr (b) {
      std::array<T, Dim> m(tmp);
      for (int i = 0; i < Dim; i++)
        tmp[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                    T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(tmp.data());
    }
    return std::inner_product(_x, _x + Dim, _x, T(0));
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<2, Dim, T, b, asymmetry> {
  static constexpr T L = -10;
  static constexpr T U = 10;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++)
      tmp[i] = x[i] - shift[i];
    auto _x = std::begin(tmp);
    if constexpr (b) {
      std::array<T, Dim> m(tmp);
      for (int i = 0; i < Dim; i++)
        tmp[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                    T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(tmp.data());
    }
    return std::accumulate(_x, _x + Dim, T(0),
                           [](auto l, auto r) { return l + std::abs(r); }) +
           std::abs(std::accumulate(_x, _x + Dim, T(1), std::multiplies<T>()));
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<3, Dim, T, b, asymmetry> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    T tmp[Dim];
    std::partial_sum(_x, _x + Dim, tmp);
    return std::inner_product(tmp, tmp + Dim, tmp, T(0));
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<4, Dim, T, b, asymmetry> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    return std::abs(*std::max_element(_x, _x + Dim, [](auto l, auto r) {
      return std::abs(l) < std::abs(r);
    }));
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<5, Dim, T, b, asymmetry> {
  static constexpr T L = -30;
  static constexpr T U = 30;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    T sum(0);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    for (int i = 0; i < Dim; i++)
      t[i] += 1;
    for (int i = 0; i < Dim - 1; i++)
      sum += (t[i] - 1) * (t[i] - 1) +
             100 * (t[i + 1] - t[i] * t[i]) * (t[i + 1] - t[i] * t[i]);
    return sum;
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<6, Dim, T, b, asymmetry> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    return std::inner_product(
        _x, _x + Dim, _x, T(0), std::plus<T>(), [](auto l, auto r) {
          return std::floor(l + T(0.5)) * std::floor(r + T(0.5));
        });
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<7, Dim, T, b, asymmetry> {
  static constexpr T L = -1.28;
  static constexpr T U = 1.28;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    constexpr T k = T(1) / T(1 + 0x7fff);
    static tool::simple_rand sr(0);
    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    T sum(0);
    for (int i = 0; i < Dim; i++)
      sum += i * t[i] * t[i] * t[i] * t[i];
    return sum + k * sr();
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<8, Dim, T, b, asymmetry> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *y, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> z;
    for (int i = 0; i < Dim; i++)
      z[i] = 10 * (y[i] - shift[i]);
    if constexpr (b) {
      std::array<T, Dim> m(z);
      for (int i = 0; i < Dim; i++)
        z[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(z.data());
    }
    T sum(0);
    for (int i = 0; i < Dim; i++) {
      T x(z[i] + T(4.209687462275036e2));
      if (x > 500) {
        sum -= (500 - std::fmod(x, T(500))) *
               std::sin(std::sqrt(std::abs(500 - std::fmod(x, T(500)))));
        T tmp = (x - 500) / T(100);
        sum += tmp * tmp / Dim;
      } else if (x < -500) {
        sum -= (-500 + std::fmod(std::abs(x), 500)) *
               std::sin(std::sqrt(std::abs(500 - std::fmod(std::abs(x), 500))));
        T tmp = (x + 500) / 100;
        sum += tmp * tmp / Dim;
      } else
        sum -= x * std::sin(std::sqrt(std::abs(x)));
    }
    sum += T(4.189828872724338e2) * Dim;

    return sum;
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<9, Dim, T, b, asymmetry> {
  static constexpr T L = -5.12;
  static constexpr T U = 5.12;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    return 10 * Dim + std::accumulate(_x, _x + Dim, T(0), [](auto l, auto r) {
             return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
           });
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<10, Dim, T, b, asymmetry> {
  static constexpr T L = -32;
  static constexpr T U = 32;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    return 20 + std::numbers::e_v<T> -
           std::exp(std::accumulate(
                        _x, _x + Dim, T(0),
                        [](auto l, auto r) {
                          return l + std::cos(2 * std::numbers::pi_v<T> * r);
                        }) /
                    Dim) -
           20 * std::exp(T(-0.2) *
                         std::sqrt(std::inner_product(_x, _x + Dim, _x, T(0)) /
                                   Dim));
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<11, Dim, T, b, asymmetry> {
  static constexpr T L = -600;
  static constexpr T U = 600;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    auto v = 1 + std::inner_product(_x, _x + Dim, _x, T(0)) * T(1) / T(4000);
    T pro(1);
    for (int i = 0; i < Dim; i++)
      pro *= std::cos(t[i] / std::sqrt(T(i + 1)));
    return v - pro;
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<12, Dim, T, b, asymmetry> {

  static constexpr T L = -50;
  static constexpr T U = 50;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    constexpr T pi = std::numbers::pi_v<T>;
    constexpr T a = 10;
    constexpr T k = 100;
    auto u = [](auto x) {
      if (x > a)
        return k * (x - a) * (x - a) * (x - a) * (x - a);
      if (x < -a)
        return k * (-x - a) * (-x - a) * (-x - a) * (-x - a);
      return T(0);
    };

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    for (int i = 0; i < Dim; i++)
      t[i] -= 1;
    return std::accumulate(_x, _x + Dim, T(0),
                           [&](auto l, auto r) { return l + u(r); }) +
           pi / Dim *
               (10 * tool::Pow<2>(std::sin(pi * (1 + T(0.25) * (t[0] + 1)))) +
                std::inner_product(
                    _x, _x + Dim - 1, _x + 1, T(0), std::plus<T>(),
                    [](auto l, auto r) {
                      return (l + 1) * (l + 1) *
                             (1 + 10 * tool::Pow<2>(std::sin(
                                           pi * (1 + T(0.25) * (r + 1))))) /
                             16;
                    }) +
                (t[Dim - 1] + 1) * (t[Dim - 1] + 1) / 16);
  }
};

template <int Dim, std::floating_point T, bool b, bool asymmetry>
  requires(Dim >= 1)
struct classic_problem<13, Dim, T, b, asymmetry> {

  static constexpr T L = -50;
  static constexpr T U = 50;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    constexpr T pi = std::numbers::pi_v<T>;
    constexpr T a = 5;
    constexpr T k = 100;
    auto u = [](auto x) {
      if (x > a)
        return k * (x - a) * (x - a) * (x - a) * (x - a);
      if (x < -a)
        return k * (-x - a) * (-x - a) * (-x - a) * (-x - a);
      return T(0);
    };

    std::array<T, Dim> t;
    for (int i = 0; i < Dim; i++)
      t[i] = x[i] - shift[i];
    auto _x = std::begin(t);
    if constexpr (b) {
      std::array<T, Dim> m(t);
      for (int i = 0; i < Dim; i++)
        t[i] = std::inner_product(mat[i].begin(), mat[i].end(), std::begin(m),
                                  T(0));
    }
    if constexpr (asymmetry) {
      detail::symmetry_breaking<Dim>(t.data());
    }
    for (int i = 0; i < Dim; i++)
      t[i] += 1;
    return std::accumulate(_x, _x + Dim, T(0),
                           [&](auto l, auto r) { return l + u(r); }) +
           T(0.1) * (tool::Pow<2>(std::sin(3 * pi * t[0])) +
                     std::inner_product(
                         _x, _x + Dim - 1, _x + 1, T(0), std::plus<T>(),
                         [](auto l, auto r) {
                           return (l - 1) * (l - 1) *
                                  (1 + tool::Pow<2>(std::sin(3 * pi * r)));
                         }) +
                     (t[Dim - 1] - 1) * (t[Dim - 1] - 1) *
                         (1 + tool::Pow<2>(std::sin(2 * pi * t[Dim - 1]))));
  }
};

} // namespace sevobench::yao_func
