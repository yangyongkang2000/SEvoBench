#pragma once
#include "cec_basical_func.hpp"

namespace sevobench::ieee_cec_set {

namespace cec2020_data {

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
using Matrix = std::vector<std::array<T, Dim>>;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<T, Dim>> shift_vector_data(7);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> rotate_matrix_data(7);

template <int Dim>
  requires(Dim >= 1)
std::array<std::array<int, Dim>, 3> perm_data;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 3 * Dim> shift_o8;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat8(4);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 4 * Dim> shift_o9;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat9(4);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 5 * Dim> shift_o10;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat10(5);

} // namespace cec2020_data

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
inline void init_cec2020_problem() noexcept {
  auto rand_perm = [](auto &m) {
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::iota(m.begin(), m.end(), int(0));
    std::shuffle(m.begin(), m.end(), gen);
  };
  for (auto &_ : cec2020_data::perm_data<Dim>)
    rand_perm(_);
  auto rand_o = [](auto &m) {
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<T> o1(T(-80), T(80));
    std::generate(m.begin(), m.end(), [&] { return o1(gen); });
  };
  for (auto &_ : cec2020_data::shift_vector_data<T, Dim>)
    rand_o(_);
  rand_o(cec2020_data::shift_o8<T, Dim>);
  rand_o(cec2020_data::shift_o9<T, Dim>);
  rand_o(cec2020_data::shift_o10<T, Dim>);
  for (auto &_ : cec2020_data::rotate_matrix_data<T, Dim>)
    tool::orthogonal_basis<Dim, T>(_);
  auto f = [](auto &m) {
    for (auto &_ : m)
      tool::orthogonal_basis<Dim, T>(_);
  };
  f(cec2020_data::mat8<T, Dim>);
  f(cec2020_data::mat9<T, Dim>);
  f(cec2020_data::mat10<T, Dim>);
}

template <int Prob_Index, int Dim, std::floating_point T>
struct cec2020_problem;

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<1, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return bent_cigar_func<Dim>(y.data()) + T(100);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<2, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return schwefel_func<Dim>(y.data()) + T(1100);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<3, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    return bi_rastrigin_func<Dim>(x, shift.data(), mat) + T(700);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<4, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return grie_rosen_func<Dim>(y.data()) + T(1900);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<5, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    constexpr int N1 = 3 * Dim / 10;
    constexpr int N2 = 3 * Dim / 10;
    constexpr int N3 = Dim - N1 - N2;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    T z[Dim];
    for (int i = 0; i < Dim; i++)
      z[i] = y[cec2020_data::perm_data<Dim>[0][i]];
    return schwefel_func<N1>(z) + rastrigin_func<N2>(z + N1) +
           ellips_func<N3>(z + N1 + N2) + T(1700);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<6, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    constexpr int N1 = Dim / 5;
    constexpr int N2 = Dim / 5;
    constexpr int N3 = 3 * Dim / 10;
    constexpr int N4 = Dim - N1 - N2 - N3;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    T z[Dim];
    for (int i = 0; i < Dim; i++)
      z[i] = y[cec2020_data::perm_data<Dim>[1][i]];
    return escaffer6_func<N1>(z) + hgbat_func<N2>(z + N1) +
           rosenbrock_func<N3>(z + N1 + N2) +
           schwefel_func<N4>(z + N1 + N2 + N3) + T(1600);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<7, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    constexpr int N1 = Dim / 10;
    constexpr int N2 = Dim / 5;
    constexpr int N3 = Dim / 5;
    constexpr int N4 = Dim / 5;
    constexpr int N5 = Dim - N1 - N2 - N3 - N4;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    T z[Dim];
    for (int i = 0; i < Dim; i++)
      z[i] = y[cec2020_data::perm_data<Dim>[2][i]];
    return escaffer6_func<N1>(z) + hgbat_func<N2>(z + N1) +
           rosenbrock_func<N3>(z + N1 + N2) +
           schwefel_func<N4>(z + N1 + N2 + N3) +
           ellips_func<N5>(z + N1 + N2 + N3 + N4) + T(2100);
  }
};
template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<8, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    const auto &shift = cec2020_data::shift_o8<T, Dim>;
    const auto &mat = cec2020_data::mat8<T, Dim>;
    constexpr T lambda[3] = {1, 10, 1};
    constexpr T delta[3] = {10, 20, 30};
    T fit[3];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat[0]);
    fit[0] = lambda[0] * rastrigin_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), shift.data() + Dim, mat[1]);
    fit[1] = lambda[1] * griewank_func<Dim>(y.data()) + 100;
    sr_func<Dim>(x, y.data(), shift.data() + 2 * Dim, mat[2]);
    fit[2] = lambda[2] * schwefel_func<Dim>(y.data()) + 200;
    return cf_cal<Dim, 3>(x, shift.data(), fit, (T *)delta) + T(2200);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<9, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    const auto &shift = cec2020_data::shift_o9<T, Dim>;
    const auto &mat = cec2020_data::mat9<T, Dim>;
    constexpr T lambda[4] = {10, T(1e-6), 10, 1};
    constexpr T delta[4] = {10, 20, 30, 40};
    T fit[4];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat[0]);
    fit[0] = lambda[0] * ackley_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), shift.data() + Dim, mat[1]);
    fit[1] = lambda[1] * ellips_func<Dim>(y.data()) + 100;
    sr_func<Dim>(x, y.data(), shift.data() + 2 * Dim, mat[2]);
    fit[2] = lambda[2] * griewank_func<Dim>(y.data()) + 200;
    sr_func<Dim>(x, y.data(), shift.data() + 3 * Dim, mat[3]);
    fit[3] = lambda[3] * rastrigin_func<Dim>(y.data()) + 300;
    return cf_cal<Dim, 4>(x, shift.data(), fit, (T *)delta) + T(2400);
  }
};
template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2020_problem<10, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    const auto &shift = cec2020_data::shift_o10<T, Dim>;
    const auto &mat = cec2020_data::mat10<T, Dim>;
    constexpr T lam[5] = {T(10), T(1), T(10), T(1e-6), T(1)};
    constexpr T delta[5] = {10, 20, 30, 40, 50};
    T fit[5];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat[0]);
    fit[0] = lam[0] * rastrigin_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), shift.data() + Dim, mat[1]);
    fit[1] = lam[1] * happycat_func<Dim>(y.data()) + 100;
    sr_func<Dim>(x, y.data(), shift.data() + 2 * Dim, mat[2]);
    fit[2] = lam[2] * ackley_func<Dim>(y.data()) + 200;
    sr_func<Dim>(x, y.data(), shift.data() + 3 * Dim, mat[3]);
    fit[3] = lam[3] * discus_func<Dim>(y.data()) + 300;
    sr_func<Dim>(x, y.data(), shift.data() + 4 * Dim, mat[4]);
    fit[4] = lam[4] * rosenbrock_func<Dim>(y.data()) + 400;
    return cf_cal<Dim, 5>(x, shift.data(), fit, (T *)delta) + T(2500);
  }
};

} // namespace sevobench::ieee_cec_set
