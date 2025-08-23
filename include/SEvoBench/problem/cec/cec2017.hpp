#pragma once

#include "cec_problem.hpp"

namespace sevobench::problem {
template <int Prob_Index, int Dim, std::floating_point T> class cec2017;

template <int Dim, std::floating_point T>
class cec2017<1, Dim, T> : public cec_common<1, Dim, T, cec2017> {
public:
  using cec_common<1, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::bent_cigar_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<2, Dim, T> : public cec_common<2, Dim, T, cec2017> {
public:
  using cec_common<2, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::sum_diff_pow_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<3, Dim, T> : public cec_common<3, Dim, T, cec2017> {
public:
  using cec_common<3, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::zakharov_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<4, Dim, T> : public cec_common<4, Dim, T, cec2017> {
public:
  using cec_common<4, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::rosenbrock_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<5, Dim, T> : public cec_common<5, Dim, T, cec2017> {
public:
  using cec_common<5, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::rastrigin_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<6, Dim, T> : public cec_common<6, Dim, T, cec2017> {
public:
  using cec_common<6, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::schwefel_F7_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<7, Dim, T> : public cec_common<7, Dim, T, cec2017> {
public:
  using cec_common<7, Dim, T, cec2017>::cec_common;
  auto operator()(std::span<const T> x) const {
    return cec_detail::bi_rastrigin_func<Dim, T>(x, this->shift, this->matrix) +
           +this->optimum_num();
  }
  static constexpr auto is_basic() { return true; }
};
template <int Dim, std::floating_point T>
class cec2017<8, Dim, T> : public cec_common<8, Dim, T, cec2017> {
public:
  using cec_common<8, Dim, T, cec2017>::cec_common;
  auto operator()(std::span<const T> x) const {
    std::array<T, Dim> y;
    std::array<T, Dim> z;
    for (int i = 0; i < Dim; i++) {
      y[i] = x[i];
      if (std::abs(y[i] - this->shift[i]) > T(0.5))
        y[i] = this->shift[i] +
               std::floor(2 * (y[i] - this->shift[i]) + T(0.5)) * T(0.5);
    }
    cec_detail::sr_func<Dim, false, T>(y, z, this->shift, this->matrix);
    return cec_detail::rastrigin_func<T, Dim>(z) + this->optimum_num();
  }
  static constexpr auto is_basic() { return true; }
};
template <int Dim, std::floating_point T>
class cec2017<9, Dim, T> : public cec_common<9, Dim, T, cec2017> {
public:
  using cec_common<9, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) { return cec_detail::levy_func(x); }
};
template <int Dim, std::floating_point T>
class cec2017<10, Dim, T> : public cec_common<10, Dim, T, cec2017> {
public:
  using cec_common<10, Dim, T, cec2017>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::schwefel_func(x);
  }
};
template <int Dim, std::floating_point T>
class cec2017<11, Dim, T> : public cec_common<11, Dim, T, cec2017> {
public:
  using cec_common<11, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x, std::to_array<std::pair<int, int>>({{1, 5}, {2, 5}, {2, 5}}),
        std::to_array<T (*)(std::span<T>)>({cec_detail::zakharov_func<T>,
                                            cec_detail::rosenbrock_func<T>,
                                            cec_detail::rastrigin_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2017<12, Dim, T> : public cec_common<12, Dim, T, cec2017> {
public:
  using cec_common<12, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x, std::to_array<std::pair<int, int>>({{3, 10}, {3, 10}, {2, 5}}),
        std::to_array<T (*)(std::span<T>)>({cec_detail::ellips_func<T>,
                                            cec_detail::schwefel_func<T>,
                                            cec_detail::bent_cigar_func<T>}));
  }
};

template <int Dim, std::floating_point T>
class cec2017<13, Dim, T> : public cec_common<13, Dim, T, cec2017> {
public:
  using cec_common<13, Dim, T, cec2017>::cec_common;
  auto operator()(std::span<const T> x) const {
    std::array<T, Dim> y;
    cec_detail::sr_func<Dim, false, T>(x, y, this->shift, this->matrix);
    std::array<T, Dim> z;
    for (int i = 0; i < Dim; i++)
      z[i] = y[this->shuffle[i]];
    for (int i = 3 * Dim / 5; i < Dim; i++) {
      z[i] = T(0.2) * T(this->shift[this->shuffle[i]] > 0 ? 1 : -1) * z[i];
    }
    return cec_detail::calculate_hybrid<Dim, T>(
               z,
               std::to_array<std::pair<int, int>>({{3, 10}, {3, 10}, {2, 5}}),
               std::to_array<T (*)(std::span<T>)>(
                   {cec_detail::bent_cigar_func<T>,
                    cec_detail::rosenbrock_func<T>,
                    cec_detail::bi_rastrigin_func<T>})) +
           T(1300);
  }
  static constexpr auto is_hybrid() { return true; }
};

template <int Dim, std::floating_point T>
class cec2017<14, Dim, T> : public cec_common<14, Dim, T, cec2017> {
public:
  using cec_common<14, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x, std::to_array<std::pair<int, int>>({{1, 5}, {1, 5}, {1, 5}, {2, 5}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::ellips_func<T>, cec_detail::ackley_func<T>,
             cec_detail::schwefel_F7_func<T>, cec_detail::rastrigin_func<T>}));
  }
};

template <int Dim, std::floating_point T>
class cec2017<15, Dim, T> : public cec_common<15, Dim, T, cec2017> {
public:
  using cec_common<15, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>({{1, 5}, {1, 5}, {3, 10}, {3, 10}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::bent_cigar_func<T>, cec_detail::hgbat_func<T>,
             cec_detail::rastrigin_func<T>, cec_detail::rosenbrock_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2017<16, Dim, T> : public cec_common<16, Dim, T, cec2017> {
public:
  using cec_common<16, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>({{1, 5}, {1, 5}, {3, 10}, {3, 10}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::escaffer6_func<T>, cec_detail::hgbat_func<T>,
             cec_detail::rosenbrock_func<T>, cec_detail::schwefel_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2017<17, Dim, T> : public cec_common<17, Dim, T, cec2017> {
public:
  using cec_common<17, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 10}, {1, 5}, {1, 5}, {1, 5}, {3, 10}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::katsuura_func<T>, cec_detail::ackley_func<T>,
             cec_detail::grie_rosen_func<T>, cec_detail::schwefel_func<T>,
             cec_detail::rastrigin_func<T>}));
  }
};

template <int Dim, std::floating_point T>
class cec2017<18, Dim, T> : public cec_common<18, Dim, T, cec2017> {
public:
  using cec_common<18, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 5}, {1, 5}, {1, 5}, {1, 5}, {1, 5}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::ellips_func<T>, cec_detail::ackley_func<T>,
             cec_detail::rastrigin_func<T>, cec_detail::hgbat_func<T>,
             cec_detail::discus_func<T>}));
  }
};

template <int Dim, std::floating_point T>
class cec2017<19, Dim, T> : public cec_common<19, Dim, T, cec2017> {
public:
  using cec_common<19, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 5}, {1, 5}, {1, 5}, {1, 5}, {1, 5}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::bent_cigar_func<T>, cec_detail::rastrigin_func<T>,
             cec_detail::grie_rosen_func<T>, cec_detail::weierstrass_func<T>,
             cec_detail::escaffer6_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2017<20, Dim, T> : public cec_common<20, Dim, T, cec2017> {
public:
  using cec_common<20, Dim, T, cec2017>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 10}, {1, 10}, {1, 5}, {1, 5}, {1, 5}, {1, 5}}),
        std::to_array<T (*)(std::span<T>)>(
            {cec_detail::hgbat_func<T>, cec_detail::katsuura_func<T>,
             cec_detail::ackley_func<T>, cec_detail::rastrigin_func<T>,
             cec_detail::schwefel_func<T>, cec_detail::schwefel_F7_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2017<21, Dim, T> : public cec_common<21, Dim, T, cec2017> {
public:
  using cec_common<21, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30}), std::to_array<T>({1, 1e-6, 1}),
        std::to_array<T>({0, 100, 200}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::rosenbrock_func<T, Dim>,
             cec_detail::ellips_func<T, Dim>,
             cec_detail::rastrigin_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
};
template <int Dim, std::floating_point T>
class cec2017<22, Dim, T> : public cec_common<22, Dim, T, cec2017> {
public:
  using cec_common<22, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30}), std::to_array<T>({1, 10, 1}),
        std::to_array<T>({0, 100, 200}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::rastrigin_func<T, Dim>,
             cec_detail::griewank_func<T, Dim>,
             cec_detail::schwefel_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
};
template <int Dim, std::floating_point T>
class cec2017<23, Dim, T> : public cec_common<23, Dim, T, cec2017> {
public:
  using cec_common<23, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40}), std::to_array<T>({1, 10, 1, 1}),
        std::to_array<T>({0, 100, 200, 300}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::rosenbrock_func<T, Dim>,
             cec_detail::ackley_func<T, Dim>, cec_detail::schwefel_func<T, Dim>,
             cec_detail::rastrigin_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 4; }
};
template <int Dim, std::floating_point T>
class cec2017<24, Dim, T> : public cec_common<24, Dim, T, cec2017> {
public:
  using cec_common<24, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40}),
        std::to_array<T>({10, 1e-6, 10, 1}),
        std::to_array<T>({0, 100, 200, 300}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::ackley_func<T, Dim>, cec_detail::ellips_func<T, Dim>,
             cec_detail::griewank_func<T, Dim>,
             cec_detail::rastrigin_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 4; }
};
template <int Dim, std::floating_point T>
class cec2017<25, Dim, T> : public cec_common<25, Dim, T, cec2017> {
public:
  using cec_common<25, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40, 50}),
        std::to_array<T>({10, 1, 10, 1e-6, 1}),
        std::to_array<T>({0, 100, 200, 300, 400}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::rastrigin_func<T, Dim>,
             cec_detail::happycat_func<T, Dim>, cec_detail::ackley_func<T, Dim>,
             cec_detail::discus_func<T, Dim>,
             cec_detail::rosenbrock_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 5; }
};
template <int Dim, std::floating_point T>
class cec2017<26, Dim, T> : public cec_common<26, Dim, T, cec2017> {
public:
  using cec_common<26, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 20, 30, 40}),
        std::to_array<T>({5e-4, 1, 10, 1, 10}),
        std::to_array<T>({0, 100, 200, 300, 400}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::escaffer6_func<T, Dim>,
             cec_detail::schwefel_func<T, Dim>,
             cec_detail::griewank_func<T, Dim>,
             cec_detail::rosenbrock_func<T, Dim>,
             cec_detail::rastrigin_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 5; }
};
template <int Dim, std::floating_point T>
class cec2017<27, Dim, T> : public cec_common<27, Dim, T, cec2017> {
public:
  using cec_common<27, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40, 50, 60}),
        std::to_array<T>({10, 10, 2.5, 1e-26, 1e-6, 5e-4}),
        std::to_array<T>({0, 100, 200, 300, 400, 500}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::hgbat_func<T, Dim>, cec_detail::rastrigin_func<T, Dim>,
             cec_detail::schwefel_func<T, Dim>,
             cec_detail::bent_cigar_func<T, Dim>,
             cec_detail::ellips_func<T, Dim>,
             cec_detail::escaffer6_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 6; }
};
template <int Dim, std::floating_point T>
class cec2017<28, Dim, T> : public cec_common<28, Dim, T, cec2017> {
public:
  using cec_common<28, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40, 50, 60}),
        std::to_array<T>({10, 10, 1e-6, 1, 1, 5e-4}),
        std::to_array<T>({0, 100, 200, 300, 400, 500}),
        std::to_array<T (*)(std::span<T, Dim>)>(
            {cec_detail::ackley_func<T, Dim>, cec_detail::griewank_func<T, Dim>,
             cec_detail::discus_func<T, Dim>,
             cec_detail::rosenbrock_func<T, Dim>,
             cec_detail::happycat_func<T, Dim>,
             cec_detail::escaffer6_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 6; }
};
template <int Dim, std::floating_point T>
class cec2017<29, Dim, T> : public cec_common<29, Dim, T, cec2017> {
public:
  using cec_common<29, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, true>(
        x, std::to_array<T>({10, 30, 50}), std::to_array<T>({1, 1, 1}),
        std::to_array<T>({0, 100, 200}),
        std::to_array({cec2017<15, Dim, T>::hybrid_evaluate,
                       cec2017<16, Dim, T>::hybrid_evaluate,
                       cec2017<17, Dim, T>::hybrid_evaluate}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
  static constexpr auto is_hybrid_composition() { return true; }
};
template <int Dim, std::floating_point T>
class cec2017<30, Dim, T> : public cec_common<30, Dim, T, cec2017> {
public:
  using cec_common<30, Dim, T, cec2017>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, true>(
        x, std::to_array<T>({10, 30, 50}), std::to_array<T>({1, 1, 1}),
        std::to_array<T>({0, 100, 200}),
        std::to_array({cec2017<15, Dim, T>::hybrid_evaluate,
                       cec2017<18, Dim, T>::hybrid_evaluate,
                       cec2017<19, Dim, T>::hybrid_evaluate}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
  static constexpr auto is_hybrid_composition() { return true; }
};
} // namespace sevobench::problem
