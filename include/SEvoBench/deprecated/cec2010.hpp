#pragma once
#include "cec_basical_func_simd.hpp"
namespace sevobench::ieee_cec_set {
namespace cec2010_data {

template <std::floating_point T, int Dim>
  requires simd_type_detail::simd_dim<T, Dim>
std::vector<std::array<T, Dim>> shift_vector_data(20);

template <int Dim>
  requires(Dim >= 1)
std::array<std::array<int, Dim>, 3> perm_data;

template <std::floating_point T, int Dim>
  requires(Dim >= 1)
std::vector<std::array<std::array<T, Dim>, Dim>> rotate_matrix_data(8);

template <std::floating_point T, int Dim, int N>
std::vector<std::array<std::array<std::array<T, Dim>, Dim>, N>>
    group_rotate_matrix_data(8);

} // namespace cec2010_data
template <int Prob_Index, int Dim, std::floating_point T>
struct cec2010_problem;

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim>
struct cec2010_problem<1, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    std::array<T, Dim> tmp;
    CEC_Base_Func_Simd::shift_func<Dim>(
        x, tmp.data(), cec2010_data::shift_vector_data<T, Dim>[0].data());
    return CEC_Base_Func_Simd::elliptic<Dim>(tmp.data());
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim>
struct cec2010_problem<2, Dim, T> {
  static constexpr T L = -5;
  static constexpr T U = 5;
  T operator()(const T *x) noexcept {
    std::array<T, Dim> tmp;
    CEC_Base_Func_Simd::shift_func<Dim>(
        x, tmp.data(), cec2010_data::shift_vector_data<T, Dim>[1].data());
    return CEC_Base_Func_Simd::rastrigin<Dim>(tmp.data());
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim>
struct cec2010_problem<3, Dim, T> {
  static constexpr T L = -32;
  static constexpr T U = 32;
  T operator()(const T *x) noexcept {
    std::array<T, Dim> tmp;
    CEC_Base_Func_Simd::shift_func<Dim>(
        x, tmp.data(), cec2010_data::shift_vector_data<T, Dim>[2].data());
    return CEC_Base_Func_Simd::ackley<Dim>(tmp.data());
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<4, Dim, T> {
  static constexpr T L = -32;
  static constexpr T U = 32;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[0][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[3][index];
    }
    CEC_Base_Func_Simd::rotate_func<M>(
        tmp.data(), cec2010_data::rotate_matrix_data<T, M>[0]);
    return 1000000 * CEC_Base_Func_Simd::elliptic<M>(tmp.data()) +
           CEC_Base_Func_Simd::elliptic<Dim - M>(tmp.data() + M);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<5, Dim, T> {
  static constexpr T L = -5;
  static constexpr T U = 5;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[1][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[4][index];
    }
    CEC_Base_Func_Simd::rotate_func<M>(
        tmp.data(), cec2010_data::rotate_matrix_data<T, M>[1]);
    return 1000000 * CEC_Base_Func_Simd::rastrigin<M>(tmp.data()) +
           CEC_Base_Func_Simd::rastrigin<Dim - M>(tmp.data() + M);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<6, Dim, T> {
  static constexpr T L = -32;
  static constexpr T U = 32;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[2][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[5][index];
    }
    CEC_Base_Func_Simd::rotate_func<M>(
        tmp.data(), cec2010_data::rotate_matrix_data<T, M>[2]);
    return 1000000 * CEC_Base_Func_Simd::ackley<M>(tmp.data()) +
           CEC_Base_Func_Simd::ackley<Dim - M>(tmp.data() + M);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<7, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[3][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[6][index];
    }
    CEC_Base_Func_Simd::rotate_func<M>(
        tmp.data(), cec2010_data::rotate_matrix_data<T, M>[3]);
    return 1000000 * CEC_Base_Func_Simd::schwefel_1<M>(tmp.data()) +
           CEC_Base_Func_Simd::sphere<Dim - M>(tmp.data() + M);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<8, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[4][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[7][index];
    }
    return 1000000 * CEC_Base_Func_Simd::rosenbrock<M>(tmp.data()) +
           CEC_Base_Func_Simd::sphere<Dim - M>(tmp.data() + M);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<9, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    constexpr int N = Dim / (2 * M);
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[5][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[8][index];
    }
    T sum(0);
    auto &mats = cec2010_data::group_rotate_matrix_data<T, M, N>[0];
    for (int i = 0; i < N; i++) {
      CEC_Base_Func_Simd::rotate_func<M>(tmp.data() + i * M, mats[i]);
      sum += CEC_Base_Func_Simd::elliptic<M>(tmp.data() + i * M);
    }
    return sum + CEC_Base_Func_Simd::elliptic<Dim / 2>(tmp.data() + Dim / 2);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<10, Dim, T> {
  static constexpr T L = -5;
  static constexpr T U = 5;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    constexpr int N = Dim / (2 * M);
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[6][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[9][index];
    }
    T sum(0);
    auto &mats = cec2010_data::group_rotate_matrix_data<T, M, N>[1];
    for (int i = 0; i < N; i++) {
      CEC_Base_Func_Simd::rotate_func<M>(tmp.data() + i * M, mats[i]);
      sum += CEC_Base_Func_Simd::rastrigin<M>(tmp.data() + i * M);
    }
    return sum + CEC_Base_Func_Simd::rastrigin<Dim / 2>(tmp.data() + Dim / 2);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<11, Dim, T> {
  static constexpr T L = -32;
  static constexpr T U = 32;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    constexpr int N = Dim / (2 * M);
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[7][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[10][index];
    }
    T sum(0);
    auto &mats = cec2010_data::group_rotate_matrix_data<T, M, N>[2];
    for (int i = 0; i < N; i++) {
      CEC_Base_Func_Simd::rotate_func<M>(tmp.data() + i * M, mats[i]);
      sum += CEC_Base_Func_Simd::ackley<M>(tmp.data() + i * M);
    }
    return sum + CEC_Base_Func_Simd::ackley<Dim / 2>(tmp.data() + Dim / 2);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<12, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    constexpr int N = Dim / (2 * M);
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[8][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[11][index];
    }
    T sum(0);
    for (int i = 0; i < N; i++) {
      sum += CEC_Base_Func_Simd::schwefel_1<M>(tmp.data() + i * M);
    }
    return sum + CEC_Base_Func_Simd::sphere<Dim / 2>(tmp.data() + Dim / 2);
  }
};

template <int Dim, std::floating_point T>
  requires simd_type_detail::simd_dim<T, Dim> &&
           simd_type_detail::simd_dim<T, (Dim / 16)>
struct cec2010_problem<13, Dim, T> {
  static constexpr T L = -100;
  static constexpr T U = 100;
  T operator()(const T *x) noexcept {
    constexpr int M = Dim / 16;
    constexpr int N = Dim / (2 * M);
    std::array<T, Dim> tmp;
    for (int i = 0; i < Dim; i++) {
      auto index = cec2010_data::perm_data<Dim>[9][i];
      tmp[i] = x[index] - cec2010_data::shift_vector_data<T, Dim>[12][index];
    }
    T sum(0);
    for (int i = 0; i < N; i++) {
      sum += CEC_Base_Func_Simd::rosenbrock<M>(tmp.data() + i * M);
    }
    return sum + CEC_Base_Func_Simd::sphere<Dim / 2>(tmp.data() + Dim / 2);
  }
};

} // namespace sevobench::ieee_cec_set
