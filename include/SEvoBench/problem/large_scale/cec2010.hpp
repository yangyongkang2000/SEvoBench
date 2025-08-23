#pragma once
#include "../../common/tool.hpp"
#include "../problem.hpp"
#include "cec_base_problem_simd.hpp"
#include <vector>
namespace sevobench::problem::simd {
template <bool Using_SIMD, int W, int G, int Index, int Dim,
          std::floating_point T,
          template <bool, int, int, int, int, std::floating_point> class Drived>
  requires(Dim % G == 0) && (Dim % (2 * G) == 0)
class cec2010_common : public problem_common<Index, Dim, T> {
protected:
  std::vector<T> shift;
  std::vector<T> matrix;
  std::vector<int> shuffle;
  int ins = 0;
  using Dri = Drived<Using_SIMD, W, G, Index, Dim, T>;

public:
  cec2010_common(int _ins) : ins(_ins) {
    constexpr auto id = Dri::problem_type_id();
    constexpr auto lb = Dri::lower_bound();
    constexpr auto ub = Dri::upper_bound();
    std::random_device rd{};
    std::default_random_engine gen{rd()};
    std::uniform_real_distribution<T> dis(lb, ub);
    shift.resize(Dim);
    std::generate_n(shift.data(), Dim, [&] { return dis(gen); });
    if constexpr (id != 1 && id != 5) {
      shuffle.resize(Dim);
      std::iota(shuffle.begin(), shuffle.end(), int(0));
      std::shuffle(shuffle.begin(), shuffle.end(), gen);
      if constexpr (need_rotate()) {
        matrix.resize(total() * G * G);
        for (int i = 0; i < total(); i++) {
          tool::generate_rotate_vector<G>(
              std::span<T>(matrix.data() + i * G * G, G * G));
        }
      }
    }
  }
  auto instance() const noexcept { return ins; }
  constexpr static auto optimum_num() { return T(0); }
  auto optimum() const noexcept {
    solution<T> s(shift.begin(), shift.end());
    s.set_fitness(optimum_num());
    return s;
  }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = Index,
                           .instance = ins,
                           .dim = Dim,
                           .lb = Dri::lower_bound(),
                           .ub = Dri::upper_bound(),
                           .optimum = optimum()};
  }

  auto operator()(std::span<const T> x) const noexcept {
    constexpr auto id = Dri::problem_type_id();
    std::array<T, Dim> tmp;
    if constexpr (id == 1 || id == 5) {
      shift_func<Dim, T>(x.data(), tmp.data(), shift.data());
    } else {
      for (int i = 0; i < Dim; i++) {
        auto idx = shuffle[i];
        tmp[idx] = x[idx] - shift[idx];
      }
      if constexpr (need_rotate()) {
        for (int i = 0; i < total(); i++)
          rotate_func<G, T>(tmp.data() + i * G, matrix.data() + i * G * G);
      }
    }
    return Dri::evaluate(tmp);
  }

  constexpr static auto need_rotate() noexcept { return true; }
  constexpr static auto problem_type_id() { return Dri::problem_type_id(); }

  constexpr static auto total() {
    constexpr auto id = problem_type_id();
    return id == 2 ? 1 : id == 3 ? (Dim / (2 * G)) : Dim / (G);
  }
};

namespace cec2010_detail {

template <int total, int G, int Dim, std::floating_point T>
inline auto group_calc(auto f, [[maybe_unused]] auto g, const T *x) noexcept {
  T sum(0);
  for (int i = 0; i < total; i++) {
    sum += f(x + i * G);
  }
  if constexpr (total * G < Dim) {
    sum += g(x + total * G);
  }
  return sum;
}

} // namespace cec2010_detail
template <bool Using_SIMD, int W, int G, int Index, int Dim,
          std::floating_point T>
class cec2010_template;
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 1, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 1, Dim, T, cec2010_template> {
public:
  using cec2010_common<Using_SIMD, W, G, 1, Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return elliptic<Using_SIMD, Dim, W, T>(x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 1; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 2, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 2, Dim, T, cec2010_template> {
public:
  using cec2010_common<Using_SIMD, W, G, 2, Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return rastrigin<Using_SIMD, Dim, W, T>(x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-5); }
  constexpr static auto upper_bound() noexcept { return T(5); }
  constexpr static auto problem_type_id() { return 1; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 3, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 3, Dim, T, cec2010_template> {
public:
  using cec2010_common<Using_SIMD, W, G, 3, Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return ackley<Using_SIMD, Dim, W, T>(x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-32); }
  constexpr static auto upper_bound() noexcept { return T(32); }
  constexpr static auto problem_type_id() { return 1; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 4, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 4, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 4, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_\
common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        [](const T *x) {
          return T(1000000) * elliptic<Using_SIMD, G, W, T>(x);
        },
        elliptic<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  static constexpr auto problem_type_id() { return 2; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 5, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 5, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 5, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        [](const T *x) {
          return T(1000000) * rastrigin<Using_SIMD, G, W, T>(x);
        },
        rastrigin<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-5); }
  constexpr static auto upper_bound() noexcept { return T(5); }
  constexpr static auto problem_type_id() { return 2; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 6, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 6, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 6, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        [](const T *x) { return T(1000000) * ackley<Using_SIMD, G, W, T>(x); },
        ackley<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-32); }
  constexpr static auto upper_bound() noexcept { return T(32); }
  constexpr static auto problem_type_id() { return 2; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 7, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 7, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 7, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        [](const T *x) {
          return T(1000000) * schwefel_1<Using_SIMD, G, W, T>(x);
        },
        sphere<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 2; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 8, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 8, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 8, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        [](const T *x) {
          return T(1000000) * rosenbrock<Using_SIMD, G, W, T>(x);
        },
        sphere<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 2; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 9, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 9, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 9, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        elliptic<Using_SIMD, G, W, T>,
        elliptic<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 3; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 10, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 10, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 10, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        rastrigin<Using_SIMD, G, W, T>,
        rastrigin<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-5); }
  constexpr static auto upper_bound() noexcept { return T(5); }
  constexpr static auto problem_type_id() { return 3; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 11, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 11, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 11, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        ackley<Using_SIMD, G, W, T>,
        ackley<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-32); }
  constexpr static auto upper_bound() noexcept { return T(32); }
  constexpr static auto problem_type_id() { return 3; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 12, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 12, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 12, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        schwefel_1<Using_SIMD, G, W, T>,
        sphere<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 3; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 13, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 13, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 13, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, 13, Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        rosenbrock<Using_SIMD, G, W, T>,
        sphere<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 3; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 14, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 14, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 14, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        elliptic<Using_SIMD, G, W, T>,
        elliptic<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 4; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 15, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 15, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 15, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        rastrigin<Using_SIMD, G, W, T>,
        rastrigin<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-5); }
  constexpr static auto upper_bound() noexcept { return T(5); }
  constexpr static auto problem_type_id() { return 4; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 16, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 16, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 16, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        ackley<Using_SIMD, G, W, T>,
        ackley<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-32); }
  constexpr static auto upper_bound() noexcept { return T(32); }
  constexpr static auto problem_type_id() { return 4; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 17, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 17, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 17, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        schwefel_1<Using_SIMD, G, W, T>,
        schwefel_1<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 4; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 18, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 18, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 18, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return cec2010_detail::group_calc<Base::total(), G, Dim, T>(
        rosenbrock<Using_SIMD, G, W, T>,
        rosenbrock<Using_SIMD, Dim - Base::total() * G, W, T>, x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 4; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 19, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 19, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 19, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return schwefel_1<Using_SIMD, Dim, W, T>(x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 5; }
  constexpr static auto need_rotate() { return false; }
};
template <bool Using_SIMD, int W, int G, int Dim, std::floating_point T>
class cec2010_template<Using_SIMD, W, G, 20, Dim, T>
    : public cec2010_common<Using_SIMD, W, G, 20, Dim, T, cec2010_template> {
  using Base = cec2010_common<Using_SIMD, W, G, 20, Dim, T, cec2010_template>;

public:
  using cec2010_common<Using_SIMD, W, G, Base::index(), Dim, T,
                       cec2010_template>::cec2010_common;
  static auto evaluate(std::span<const T, Dim> x) noexcept {
    return rosenbrock<Using_SIMD, Dim, W, T>(x.data());
  }
  constexpr static auto lower_bound() noexcept { return T(-100); }
  constexpr static auto upper_bound() noexcept { return T(100); }
  constexpr static auto problem_type_id() { return 5; }
  constexpr static auto need_rotate() { return false; }
};
template <int G, bool Using_SIMD = (simd_id() > 0), int W = 8>
struct cec2010_setting {
  template <int Index, int Dim, std::floating_point T>
  using cec2010 =
      cec2010_template<Using_SIMD, Using_SIMD ? std::min(simd_width<T>, W) : -1,
                       G, Index, Dim, T>;
};
} // namespace sevobench::problem::simd
