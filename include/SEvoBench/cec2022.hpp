#pragma once
#include "cec_basical_func.hpp"

namespace sevobench::ieee_cec_set {
namespace cec2022_data {

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<T, Dim>> shift_vector_data(8);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> rotate_matrix_data(8);

template <int Dim> std::array<std::array<int, Dim>, 3> perm_data;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 5 * Dim> shift_o9;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat9(5);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 3 * Dim> shift_o10;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat10(3);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 5 * Dim> shift_o11;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat11(5);

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::array<T, 6 * Dim> shift_o12;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> mat12(6);

} // namespace cec2022_data

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
inline void init_cec2022_problem() noexcept {
  auto rand_perm = [](auto &m) {
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::iota(m.begin(), m.end(), int(0));
    std::shuffle(m.begin(), m.end(), gen);
  };
  for (auto &_ : cec2022_data::perm_data<Dim>)
    rand_perm(_);
  auto rand_o = [](auto &m) {
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<T> o1(T(-80), T(80));
    std::generate(m.begin(), m.end(), [&] { return o1(gen); });
  };
  for (auto &_ : cec2022_data::shift_vector_data<T, Dim>)
    rand_o(_);
  rand_o(cec2022_data::shift_o9<T, Dim>);
  rand_o(cec2022_data::shift_o10<T, Dim>);
  rand_o(cec2022_data::shift_o11<T, Dim>);
  rand_o(cec2022_data::shift_o12<T, Dim>);
  for (auto &_ : cec2022_data::rotate_matrix_data<T, Dim>)
    tool::orthogonal_basis<Dim, T>(_);
  auto f = [](auto &m) {
    for (auto &_ : m)
      tool::orthogonal_basis<Dim, T>(_);
  };
  f(cec2022_data::mat9<T, Dim>);
  f(cec2022_data::mat10<T, Dim>);
  f(cec2022_data::mat11<T, Dim>);
  f(cec2022_data::mat12<T, Dim>);
}

template <int Prob_Index, int Dim, std::floating_point T>
struct cec2022_problem;

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<1, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return zakharov_func<Dim>(y.data()) + T(300);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<2, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return rosenbrock_func<Dim>(y.data()) + T(400);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<3, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return escaffer6_func<Dim>(y.data()) + T(600);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<4, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return step_rastrigin_func<Dim>(y.data()) + T(800);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<5, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    return levy_func<Dim>(y.data()) + T(900);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<6, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    constexpr int N1 = 2 * Dim / 5;
    constexpr int N2 = 2 * Dim / 5;
    constexpr int N3 = Dim - N1 - N2;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    T z[Dim];
    for (int i = 0; i < Dim; i++)
      z[i] = y[cec2022_data::perm_data<Dim>[0][i]];
    return bent_cigar_func<N1>(z) + hgbat_func<N2>(z + N1) +
           rastrigin_func<N3>(z + N1 + N2) + T(1800);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<7, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    constexpr int N1 = Dim / 10;
    constexpr int N2 = Dim / 5;
    constexpr int N3 = Dim / 5;
    constexpr int N4 = Dim / 5;
    constexpr int N5 = Dim / 10;
    constexpr int N6 = Dim - N1 - N2 - N3 - N4 - N5;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    T z[Dim];
    for (int i = 0; i < Dim; i++)
      z[i] = y[cec2022_data::perm_data<Dim>[1][i]];
    return hgbat_func<N1>(z) + katsuura_func<N2>(z + N1) +
           ackley_func<N3>(z + N1 + N2) + rastrigin_func<N4>(z + N1 + N2 + N3) +
           schwefel_func<N5>(z + N1 + N2 + N3 + N4) +
           schwefel_F7_func<N6>(z + N1 + N2 + N3 + N4 + N5) + T(2000);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<8, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x, const auto &shift, const auto &mat) noexcept {
    using namespace CEC_Base_Func;
    constexpr int N1 = 3 * Dim / 10;
    constexpr int N2 = Dim / 5;
    constexpr int N3 = Dim / 5;
    constexpr int N4 = Dim / 10;
    constexpr int N5 = Dim - N1 - N2 - N3 - N4;
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), shift.data(), mat);
    T z[Dim];
    for (int i = 0; i < Dim; i++)
      z[i] = y[cec2022_data::perm_data<Dim>[2][i]];
    return katsuura_func<N1>(z) + happycat_func<N2>(z + N1) +
           grie_rosen_func<N3>(z + N1 + N2) +
           schwefel_func<N4>(z + N1 + N2 + N3) +
           ackley_func<N5>(z + N1 + N2 + N3 + N4) + T(2200);
  }
};
template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<9, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    constexpr T lam[5] = {T(1), T(1e-6), T(1e-6), T(1e-6), T(1e-6)};
    constexpr T delta[5] = {10, 20, 30, 40, 50};
    T fit[5];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o9<T, Dim>.data(),
                 cec2022_data::mat9<T, Dim>[0]);
    fit[0] = lam[0] * rosenbrock_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o9<T, Dim>.data() + Dim,
                 cec2022_data::mat9<T, Dim>[1]);
    fit[1] = lam[1] * ellips_func<Dim>(y.data()) + 200;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o9<T, Dim>.data() + 2 * Dim,
                 cec2022_data::mat9<T, Dim>[2]);
    fit[2] = lam[2] * bent_cigar_func<Dim>(y.data()) + 300;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o9<T, Dim>.data() + 3 * Dim,
                 cec2022_data::mat9<T, Dim>[3]);
    fit[3] = lam[3] * discus_func<Dim>(y.data()) + 100;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o9<T, Dim>.data() + 4 * Dim,
                 cec2022_data::mat9<T, Dim>[4]);
    fit[4] = lam[4] * ellips_func<Dim>(y.data()) + 400;
    return cf_cal<Dim, 5>(x, cec2022_data::shift_o9<T, Dim>.data(), fit,
                          (T *)delta) +
           T(2300);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<10, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    constexpr T delta[3] = {20, 10, 10};
    T fit[3];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o10<T, Dim>.data(),
                 cec2022_data::mat10<T, Dim>[0]);
    fit[0] = schwefel_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o10<T, Dim>.data() + Dim,
                 cec2022_data::mat10<T, Dim>[1]);
    fit[1] = rastrigin_func<Dim>(y.data()) + 200;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o10<T, Dim>.data() + 2 * Dim,
                 cec2022_data::mat10<T, Dim>[2]);
    fit[2] = hgbat_func<Dim>(y.data()) + 100;
    return cf_cal<Dim, 3>(x, cec2022_data::shift_o10<T, Dim>.data(), fit,
                          (T *)delta) +
           T(2400);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<11, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    constexpr T lam[5] = {T(5e-4), T(1), T(10), T(1), T(10)};
    constexpr T delta[5] = {20, 20, 30, 30, 20};
    T fit[5];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o11<T, Dim>.data(),
                 cec2022_data::mat11<T, Dim>[0]);
    fit[0] = lam[0] * escaffer6_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o11<T, Dim>.data() + Dim,
                 cec2022_data::mat11<T, Dim>[1]);
    fit[1] = lam[1] * schwefel_func<Dim>(y.data()) + 200;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o11<T, Dim>.data() + 2 * Dim,
                 cec2022_data::mat11<T, Dim>[2]);
    fit[2] = lam[2] * griewank_func<Dim>(y.data()) + 300;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o11<T, Dim>.data() + 3 * Dim,
                 cec2022_data::mat11<T, Dim>[3]);
    fit[3] = lam[3] * rosenbrock_func<Dim>(y.data()) + 400;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o11<T, Dim>.data() + 4 * Dim,
                 cec2022_data::mat11<T, Dim>[4]);
    fit[4] = lam[4] * rastrigin_func<Dim>(y.data()) + 200;
    return cf_cal<Dim, 5>(x, cec2022_data::shift_o11<T, Dim>.data(), fit,
                          (T *)delta) +
           T(2600);
  }
};

template <int Dim, std::floating_point T>
  requires(Dim >= 1)
struct cec2022_problem<12, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    using namespace CEC_Base_Func;
    constexpr T lam[6] = {T(10), T(10), T(2.5), T(1e-26), T(1e-6), T(5e-4)};
    constexpr T delta[6] = {10, 20, 30, 40, 50, 60};
    T fit[6];
    std::array<T, Dim> y;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o12<T, Dim>.data(),
                 cec2022_data::mat12<T, Dim>[0]);
    fit[0] = lam[0] * hgbat_func<Dim>(y.data());
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o12<T, Dim>.data() + Dim,
                 cec2022_data::mat12<T, Dim>[1]);
    fit[1] = lam[1] * rastrigin_func<Dim>(y.data()) + 300;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o12<T, Dim>.data() + 2 * Dim,
                 cec2022_data::mat12<T, Dim>[2]);
    fit[2] = lam[2] * schwefel_func<Dim>(y.data()) + 500;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o12<T, Dim>.data() + 3 * Dim,
                 cec2022_data::mat12<T, Dim>[3]);
    fit[3] = lam[3] * bent_cigar_func<Dim>(y.data()) + 100;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o12<T, Dim>.data() + 4 * Dim,
                 cec2022_data::mat12<T, Dim>[4]);
    fit[4] = lam[4] * ellips_func<Dim>(y.data()) + 400;
    sr_func<Dim>(x, y.data(), cec2022_data::shift_o12<T, Dim>.data() + 5 * Dim,
                 cec2022_data::mat12<T, Dim>[5]);
    fit[5] = lam[5] * escaffer6_func<Dim>(y.data()) + 200;
    return cf_cal<Dim, 6>(x, cec2022_data::shift_o12<T, Dim>.data(), fit,
                          (T *)delta) +
           T(2700);
  }
};
} // namespace sevobench::ieee_cec_set
