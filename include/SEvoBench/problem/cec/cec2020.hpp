#pragma once
#include "cec2014.hpp"
#include "cec2017.hpp"
namespace sevobench::problem {
template <int Prob_Index, int Dim, std::floating_point T> struct cec2020;

template <int Dim, std::floating_point T>
struct cec2020<1, Dim, T> : public cec2017<1, Dim, T> {
  using cec2017<1, Dim, T>::cec2017;

  cec2020(const std::string &dir_name)
      : cec2017<1, Dim, T>(dir_name, index()) {}

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
struct cec2020<2, Dim, T> : public cec2014<11, Dim, T> {
  using cec2014<11, Dim, T>::cec2014;

  cec2020(const std::string &dir_name)
      : cec2014<11, Dim, T>(dir_name, index()) {}

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
struct cec2020<3, Dim, T> : public cec2017<7, Dim, T> {
  using cec2017<7, Dim, T>::cec2017;

  cec2020(const std::string &dir_name)
      : cec2017<7, Dim, T>(dir_name, index()) {}

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
struct cec2020<4, Dim, T> : public cec_common<4, Dim, T, cec2020> {
  using cec_common<4, Dim, T, cec2020>::cec_common;
  cec2020(const std::string &dir_name)
      : cec_common<4, Dim, T, cec2020>(dir_name, 7) {}
  static constexpr auto optimum_num() { return 1900; }
  static auto evaluate(std::span<T, Dim> x) {
    return cec_detail::grie_rosen_func<T, Dim>(x);
  }
};
template <int Dim, std::floating_point T>
struct cec2020<5, Dim, T> : public cec2014<17, Dim, T> {
  using cec2014<17, Dim, T>::cec2014;

  cec2020(const std::string &dir_name) : cec2014<17, Dim, T>(dir_name, 4) {}

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
struct cec2020<6, Dim, T> : public cec2017<16, Dim, T> {
  using cec2017<16, Dim, T>::cec2017;

  cec2020(const std::string &dir_name) : cec2017<16, Dim, T>(dir_name, 16) {}

  static constexpr auto index() { return 6; }

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
struct cec2020<7, Dim, T> : public cec2014<21, Dim, T> {
  using cec2014<21, Dim, T>::cec2014;

  cec2020(const std::string &dir_name) : cec2014<21, Dim, T>(dir_name, 6) {}

  static constexpr auto index() { return 7; }

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
struct cec2020<8, Dim, T> : public cec2017<22, Dim, T> {
  using cec2017<22, Dim, T>::cec2017;

  cec2020(const std::string &dir_name) : cec2017<22, Dim, T>(dir_name, 22) {}

  static constexpr auto index() { return 8; }

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
struct cec2020<9, Dim, T> : public cec2017<24, Dim, T> {
  using cec2017<24, Dim, T>::cec2017;

  cec2020(const std::string &dir_name) : cec2017<24, Dim, T>(dir_name, 24) {}

  static constexpr auto index() { return 9; }

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
struct cec2020<10, Dim, T> : public cec2017<25, Dim, T> {
  using cec2017<25, Dim, T>::cec2017;

  cec2020(const std::string &dir_name) : cec2017<25, Dim, T>(dir_name, 25) {}

  static constexpr auto index() { return 10; }

  auto problem_information() const noexcept {
    return problem_info<T>{.index = index(),
                           .instance = this->ins,
                           .dim = this->dim(),
                           .lb = this->lower_bound(),
                           .ub = this->upper_bound(),
                           .optimum = this->optimum_solution()};
  }
};
} // namespace sevobench::problem