#pragma once
#include "cec_problem.hpp"
namespace sevobench::problem {
template <int Prob_Index, int Dim, std::floating_point T> class cec2014;
template <int Dim, std::floating_point T>
class cec2014<1, Dim, T> : public cec_common<1, Dim, T, cec2014> {
public:
  using cec_common<1, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::ellips_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<2, Dim, T> : public cec_common<2, Dim, T, cec2014> {
public:
  using cec_common<2, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::bent_cigar_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<3, Dim, T> : public cec_common<3, Dim, T, cec2014> {
public:
  using cec_common<3, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::discus_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<4, Dim, T> : public cec_common<4, Dim, T, cec2014> {
public:
  using cec_common<4, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::rosenbrock_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<5, Dim, T> : public cec_common<5, Dim, T, cec2014> {
public:
  using cec_common<5, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::ackley_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<6, Dim, T> : public cec_common<6, Dim, T, cec2014> {
public:
  using cec_common<6, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::weierstrass_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<7, Dim, T> : public cec_common<7, Dim, T, cec2014> {
public:
  using cec_common<7, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::griewank_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<8, Dim, T> : public cec_common<8, Dim, T, cec2014> {
public:
  using cec_common<8, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::rastrigin_func<T, Dim>(x);
  }
  static constexpr auto is_only_shift() { return true; }
};
template <int Dim, std::floating_point T>
class cec2014<9, Dim, T> : public cec_common<9, Dim, T, cec2014> {
public:
  using cec_common<9, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::rastrigin_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<10, Dim, T> : public cec_common<10, Dim, T, cec2014> {
public:
  using cec_common<10, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::schwefel_func<T, Dim>(x);
  }
  static constexpr auto is_only_shift() { return true; }
};
template <int Dim, std::floating_point T>
class cec2014<11, Dim, T> : public cec_common<11, Dim, T, cec2014> {
public:
  using cec_common<11, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::schwefel_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<12, Dim, T> : public cec_common<12, Dim, T, cec2014> {
public:
  using cec_common<12, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::katsuura_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<13, Dim, T> : public cec_common<13, Dim, T, cec2014> {
public:
  using cec_common<13, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::happycat_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<14, Dim, T> : public cec_common<14, Dim, T, cec2014> {
public:
  using cec_common<14, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::hgbat_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<15, Dim, T> : public cec_common<15, Dim, T, cec2014> {
public:
  using cec_common<15, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::grie_rosen_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<16, Dim, T> : public cec_common<16, Dim, T, cec2014> {
public:
  using cec_common<16, Dim, T, cec2014>::cec_common;
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::escaffer6_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
class cec2014<17, Dim, T> : public cec_common<17, Dim, T, cec2014> {
public:
  using cec_common<17, Dim, T, cec2014>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x, std::to_array<std::pair<int, int>>({{3, 10}, {3, 10}, {2, 5}}),
        std::to_array({cec_detail::schwefel_func<T>,
                       cec_detail::rastrigin_func<T>,
                       cec_detail::ellips_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2014<18, Dim, T> : public cec_common<18, Dim, T, cec2014> {
public:
  using cec_common<18, Dim, T, cec2014>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x, std::to_array<std::pair<int, int>>({{3, 10}, {3, 10}, {2, 5}}),
        std::to_array({cec_detail::bent_cigar_func<T>,
                       cec_detail::hgbat_func<T>,
                       cec_detail::rastrigin_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2014<19, Dim, T> : public cec_common<19, Dim, T, cec2014> {
public:
  using cec_common<19, Dim, T, cec2014>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>({{1, 5}, {1, 5}, {3, 10}, {3, 10}}),
        std::to_array(
            {cec_detail::griewank_func<T>, cec_detail::weierstrass_func<T>,
             cec_detail::rosenbrock_func<T>, cec_detail::escaffer6_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2014<20, Dim, T> : public cec_common<20, Dim, T, cec2014> {
public:
  using cec_common<20, Dim, T, cec2014>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>({{1, 5}, {1, 5}, {3, 10}, {3, 10}}),
        std::to_array({cec_detail::hgbat_func<T>, cec_detail::discus_func<T>,
                       cec_detail::grie_rosen_func<T>,
                       cec_detail::rastrigin_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2014<21, Dim, T> : public cec_common<21, Dim, T, cec2014> {
public:
  using cec_common<21, Dim, T, cec2014>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 10}, {1, 5}, {1, 5}, {1, 5}, {3, 10}}),
        std::to_array({cec_detail::escaffer6_func<T>, cec_detail::hgbat_func<T>,
                       cec_detail::rosenbrock_func<T>,
                       cec_detail::schwefel_func<T>,
                       cec_detail::ellips_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2014<22, Dim, T> : public cec_common<22, Dim, T, cec2014> {
public:
  using cec_common<22, Dim, T, cec2014>::cec_common;
  static auto hybrid_evaluate(std::span<T, Dim> x) {
    return cec_detail::calculate_hybrid<Dim>(
        x,
        std::to_array<std::pair<int, int>>(
            {{1, 10}, {1, 5}, {1, 5}, {1, 5}, {3, 10}}),
        std::to_array(
            {cec_detail::katsuura_func<T>, cec_detail::happycat_func<T>,
             cec_detail::grie_rosen_func<T>, cec_detail::schwefel_func<T>,
             cec_detail::ackley_func<T>}));
  }
};
template <int Dim, std::floating_point T>
class cec2014<23, Dim, T> : public cec_common<23, Dim, T, cec2014> {
public:
  using cec_common<23, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z, std::span<const int>) {
    constexpr auto lams = std::to_array<T>({1, 1e-6, 1e-26, 1e-6, 1e-6});
    constexpr auto biases = std::to_array<T>({0, 100, 200, 300, 400});
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
};
template <int Dim, std::floating_point T>
class cec2014<24, Dim, T> : public cec_common<24, Dim, T, cec2014> {
public:
  using cec_common<24, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z, std::span<const int>) {
    constexpr auto deltas = std::to_array<T>({20, 20, 20});
    constexpr auto lams = std::to_array<T>({1, 1, 1});
    constexpr auto biases = std::to_array<T>({0, 100, 200});
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
    ;
  }
  static constexpr auto cf_num() { return 3; }
};
template <int Dim, std::floating_point T>
class cec2014<25, Dim, T> : public cec_common<25, Dim, T, cec2014> {
public:
  using cec_common<25, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 30, 50}), std::to_array<T>({0.25, 1, 1e-7}),
        std::to_array<T>({0, 100, 200}),
        std::to_array({cec_detail::schwefel_func<T, Dim>,
                       cec_detail::rastrigin_func<T, Dim>,
                       cec_detail::ellips_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
};
template <int Dim, std::floating_point T>
class cec2014<26, Dim, T> : public cec_common<26, Dim, T, cec2014> {
public:
  using cec_common<26, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 10, 10, 10, 10}),
        std::to_array<T>({0.25, 1, 1e-7, 2.5, 10}),
        std::to_array<T>({0, 100, 200, 300, 400}),
        std::to_array({cec_detail::schwefel_func<T, Dim>,
                       cec_detail::happycat_func<T, Dim>,
                       cec_detail::ellips_func<T, Dim>,
                       cec_detail::weierstrass_func<T, Dim>,
                       cec_detail::griewank_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 5; }
};
template <int Dim, std::floating_point T>
class cec2014<27, Dim, T> : public cec_common<27, Dim, T, cec2014> {
public:
  using cec_common<27, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 10, 10, 20, 20}),
        std::to_array<T>({10, 10, 2.5, 25, 1e-6}),
        std::to_array<T>({0, 100, 200, 300, 400}),
        std::to_array({cec_detail::hgbat_func<T, Dim>,
                       cec_detail::rastrigin_func<T, Dim>,
                       cec_detail::schwefel_func<T, Dim>,
                       cec_detail::weierstrass_func<T, Dim>,
                       cec_detail::ellips_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 5; }
};
template <int Dim, std::floating_point T>
class cec2014<28, Dim, T> : public cec_common<28, Dim, T, cec2014> {
public:
  using cec_common<28, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, false>(
        x, std::to_array<T>({10, 20, 30, 40, 50}),
        std::to_array<T>({2.5, 10, 2.5, 5e-4, 1e-6}),
        std::to_array<T>({0, 100, 200, 300, 400}),
        std::to_array({cec_detail::grie_rosen_func<T, Dim>,
                       cec_detail::happycat_func<T, Dim>,
                       cec_detail::schwefel_func<T, Dim>,
                       cec_detail::escaffer6_func<T, Dim>,
                       cec_detail::ellips_func<T, Dim>}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 5; }
};
template <int Dim, std::floating_point T>
class cec2014<29, Dim, T> : public cec_common<29, Dim, T, cec2014> {
public:
  using cec_common<29, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, true>(
        x, std::to_array<T>({10, 30, 50}), std::to_array<T>({1, 1, 1}),
        std::to_array<T>({0, 100, 200}),
        std::to_array({cec2014<17, Dim, T>::hybrid_evaluate,
                       cec2014<18, Dim, T>::hybrid_evaluate,
                       cec2014<19, Dim, T>::hybrid_evaluate}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
  static constexpr auto is_hybrid_composition() { return true; }
};
template <int Dim, std::floating_point T>
class cec2014<30, Dim, T> : public cec_common<30, Dim, T, cec2014> {
public:
  using cec_common<30, Dim, T, cec2014>::cec_common;
  static auto composition_evaluate(std::span<const T> x, std::span<const T> y,
                                   std::span<const T> z,
                                   std::span<const int> s) {
    return cec_detail::calculate_composition<Dim, true>(
        x, std::to_array<T>({10, 30, 50}), std::to_array<T>({1, 1, 1}),
        std::to_array<T>({0, 100, 200}),
        std::to_array({cec2014<20, Dim, T>::hybrid_evaluate,
                       cec2014<21, Dim, T>::hybrid_evaluate,
                       cec2014<22, Dim, T>::hybrid_evaluate}),
        y, z, s);
  }
  static constexpr auto cf_num() { return 3; }
  static constexpr auto is_hybrid_composition() { return true; }
};

} // namespace sevobench::problem
