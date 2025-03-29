#pragma once

#include "cec2014.hpp"
#include "cec2017.hpp"

namespace sevobench::problem {
template <int Prob_Index, int Dim, std::floating_point T> struct cec2022;

template <int Dim, std::floating_point T>
struct cec2022<1, Dim, T> : public cec2017<3, Dim, T> {
  using cec2017<3, Dim, T>::cec2017;

  cec2022(const std::string &dir_name)
      : cec2017<3, Dim, T>(dir_name, index()) {}

  static constexpr auto index() { return 1; }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = index(),
                           .instance = this->ins,
                           .dim = this->dim(),
                           .lb = this->lower_bound(),
                           .ub = this->upper_bound(),
                           .optimum = this->optimum_solution()};
  }
};

template <int Dim, std::floating_point T>
struct cec2022<2, Dim, T> : public cec2017<4, Dim, T> {
  using cec2017<4, Dim, T>::cec2017;

  cec2022(const std::string &dir_name)
      : cec2017<4, Dim, T>(dir_name, index()) {}

  static constexpr auto index() { return 2; }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = index(),
                           .instance = this->ins,
                           .dim = this->dim(),
                           .lb = this->lower_bound(),
                           .ub = this->upper_bound(),
                           .optimum = this->optimum_solution()};
  }
};

template <int Dim, std::floating_point T>
struct cec2022<3, Dim, T> : public cec2017<6, Dim, T> {
  using cec2017<6, Dim, T>::cec2017;

  cec2022(const std::string &dir_name)
      : cec2017<6, Dim, T>(dir_name, index()) {}

  static constexpr auto index() { return 3; }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = index(),
                           .instance = this->ins,
                           .dim = this->dim(),
                           .lb = this->lower_bound(),
                           .ub = this->upper_bound(),
                           .optimum = this->optimum_solution()};
  }
};

template <int Dim, std::floating_point T>
struct cec2022<4, Dim, T> : public cec2017<8, Dim, T> {
  using cec2017<8, Dim, T>::cec2017;

  cec2022(const std::string &dir_name)
      : cec2017<8, Dim, T>(dir_name, index()) {}

  static constexpr auto index() { return 4; }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = index(),
                           .instance = this->ins,
                           .dim = this->dim(),
                           .lb = this->lower_bound(),
                           .ub = this->upper_bound(),
                           .optimum = this->optimum_solution()};
  }
};

template <int Dim, std::floating_point T>
struct cec2022<5, Dim, T> : public cec2017<9, Dim, T> {
  using cec2017<9, Dim, T>::cec2017;

  cec2022(const std::string &dir_name)
      : cec2017<9, Dim, T>(dir_name, index()) {}

  static constexpr auto index() { return 5; }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = index(),
                           .instance = this->ins,
                           .dim = this->dim(),
                           .lb = this->lower_bound(),
                           .ub = this->upper_bound(),
                           .optimum = this->optimum_solution()};
  }
};

template <int Dim, std::floating_point T>
struct cec2022<7, Dim, T> : public cec_common<7, Dim, T, cec2022> {
  using cec_common<7, Dim, T, cec2022>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 10}, {1, 5}, {1, 5}, {1, 5}, {1, 10}, {1, 5}}),
        std::to_array(
            {cec_detail::hgbat_func<T>, cec_detail::katsuura_func<T>,
             cec_detail::ackley_func<T>, cec_detail::rastrigin_func<T>,
             cec_detail::schwefel_func<T>, cec_detail::schwefel_F7_func<T>}));
  }
  static constexpr auto optimum_num() { return T(2000); }
};

template <int Dim, std::floating_point T>
struct cec2022<9, Dim, T> : public cec_common<9, Dim, T, cec2022> {
  using cec_common<9, Dim, T, cec2022>::cec_common;

  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z, std::span<const int>) {
    constexpr auto lams = std::to_array<T>({1, 1e-6, 1e-26, 1e-6, 1e-6});
    constexpr auto biases = std::to_array<T>({0, 200, 300, 100, 400});
    constexpr auto deltas = std::to_array<T>({10, 20, 30, 40, 50});
    constexpr auto f = std::to_array(
        {cec_detail::rosenbrock_func<T, Dim>, cec_detail::ellips_func<T, Dim>,
         cec_detail::bent_cigar_func<T, Dim>, cec_detail::discus_func<T, Dim>,
         cec_detail::ellips_func<T, Dim>});
    std::array<T, cf_num()> fits;
    std::array<T, Dim> tmp_x;
    for (int i = 0; i < cf_num(); i++) {
      if (i != cf_num() - 1) {
        cec_detail::sr_func<Dim, false, T>(x, tmp_x, y.subspan(i * Dim),
                                           z.subspan(i * Dim * Dim));
      } else {
        cec_detail::sr_func<Dim, true, T>(x, tmp_x, y.subspan(i * Dim),
                                          z.subspan(i * Dim * Dim));
      }
      fits[i] = lams[i] * f[i](tmp_x) + biases[i];
    }
    return cec_detail::cf_cal<Dim, cf_num(), T>(x, y, fits, deltas);
  }
  static constexpr auto cf_num() { return 5; }
  static constexpr auto optimum_num() { return T(2300); }
};

template <int Dim, std::floating_point T>
struct cec2022<10, Dim, T> : public cec_common<10, Dim, T, cec2022> {
  using cec_common<10, Dim, T, cec2022>::cec_common;

  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z, std::span<const int>) {
    constexpr auto deltas = std::to_array<T>({20, 10, 10});
    constexpr auto lams = std::to_array<T>({1, 1, 1});
    constexpr auto biases = std::to_array<T>({0, 200, 100});
    constexpr auto f = std::to_array({cec_detail::schwefel_func<T, Dim>,
                                      cec_detail::rastrigin_func<T, Dim>,
                                      cec_detail::hgbat_func<T, Dim>});
    std::array<T, cf_num()> fits;
    std::array<T, Dim> tmp_x;
    for (int i = 0; i < cf_num(); i++) {
      if (i != 0) {
        cec_detail::sr_func<Dim, false, T>(x, tmp_x, y.subspan(i * Dim),
                                           z.subspan(i * Dim * Dim));
      } else {
        cec_detail::sr_func<Dim, true, T>(x, tmp_x, y, z);
      }
      fits[i] = lams[i] * f[i](tmp_x) + biases[i];
    }
    return cec_detail::cf_cal<Dim, cf_num(), T>(x, y, fits, deltas);
  }
  static constexpr auto cf_num() { return 3; }
  static constexpr auto optimum_num() { return T(2400); }
};

template <int Dim, std::floating_point T>
struct cec2022<11, Dim, T> : public cec_common<11, Dim, T, cec2022> {
  using cec_common<11, Dim, T, cec2022>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({20, 20, 30, 30, 20}),
        std::to_array<T>({5e-4, 1, 10, 1, 10}),
        std::to_array<T>({0, 200, 300, 400, 200}),
        std::to_array({cec_detail::escaffer6_func<T, Dim>,
                       cec_detail::schwefel_func<T, Dim>,
                       cec_detail::griewank_func<T, Dim>,
                       cec_detail::rosenbrock_func<T, Dim>,
                       cec_detail::rastrigin_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 5; }
  static constexpr auto optimum_num() { return T(2600); }
};

template <int Dim, std::floating_point T>
struct cec2022<12, Dim, T> : public cec_common<12, Dim, T, cec2022> {
  using cec_common<12, Dim, T, cec2022>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40, 50, 60}),
        std::to_array<T>({10, 10, 2.5, 1e-26, 1e-6, 5e-4}),
        std::to_array<T>({0, 300, 500, 100, 400, 200}),
        std::to_array({cec_detail::hgbat_func<T, Dim>,
                       cec_detail::rastrigin_func<T, Dim>,
                       cec_detail::schwefel_func<T, Dim>,
                       cec_detail::bent_cigar_func<T, Dim>,
                       cec_detail::ellips_func<T, Dim>,
                       cec_detail::escaffer6_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 6; }
  static constexpr auto optimum_num() { return T(2700); }
};
template <int Dim, std::floating_point T>
struct cec2022<6, Dim, T> : public cec_common<6, Dim, T, cec2022> {
  using cec_common<6, Dim, T, cec2022>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x, std::to_array<std::pair<int, int>>({{2, 5}, {2, 5}, {1, 5}}),
        std::to_array({cec_detail::bent_cigar_func<T>,
                       cec_detail::hgbat_func<T>,
                       cec_detail::rastrigin_func<T>}));
  }
  static constexpr auto optimum_num() { return T(1800); }
};
template <int Dim, std::floating_point T>
struct cec2022<8, Dim, T> : public cec_common<8, Dim, T, cec2022> {
  using cec_common<8, Dim, T, cec2022>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{3, 10}, {1, 5}, {1, 5}, {1, 10}, {1, 5}}),
        std::to_array(
            {cec_detail::katsuura_func<T>, cec_detail::happycat_func<T>,
             cec_detail::grie_rosen_func<T>, cec_detail::schwefel_func<T>,
             cec_detail::ackley_func<T>}));
  }
  static constexpr auto optimum_num() { return T(2200); }
};
} // namespace sevobench::problem
