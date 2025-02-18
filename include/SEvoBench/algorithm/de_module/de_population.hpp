
#pragma once

#include "../../common/population.hpp"
#include "../evolutionary_algorithm.hpp"

namespace sevobench::de_module {
template <std::floating_point T> class de_population {
public:
  virtual void resize(population<T> &) {};

  virtual ~de_population() {}
};

template <std::floating_point T>
class linear_reduction final : public de_population<T> {
  const evolutionary_algorithm &alg;
  const int N_Min;
  const int N_Max;

public:
  linear_reduction(const evolutionary_algorithm &_alg, int n_min, int n_max)
      : alg(_alg), N_Min(n_min), N_Max(n_max) {}

  void resize(population<T> &pop) override {
    auto n = std::max(
        static_cast<int>(std::round(T(N_Min - N_Max) * T(alg.current_fes()) /
                                        T(alg.max_fes()) +
                                    T(N_Max))),
        N_Min);
    if (n < pop.pop_size()) {
      std::nth_element(
          pop.begin(), pop.begin() + n, pop.begin() + pop.pop_size(),
          [](auto &x, auto &y) { return x.fitness() < y.fitness(); });
      pop.resize(n);
    }
  }
};

template <std::floating_point T>
class non_linear_reduction final : public de_population<T> {
  const evolutionary_algorithm &alg;
  const int N_Min;
  const int N_Max;

public:
  non_linear_reduction(const evolutionary_algorithm &_alg, int n_min, int n_max)
      : alg(_alg), N_Min(n_min), N_Max(n_max) {}

  void resize(population<T> &pop) override {
    auto ratio = T(alg.current_fes()) / T(alg.max_fes());
    auto n = std::max(
        static_cast<int>(std::round(
            T(N_Min - N_Max) * std::pow(ratio, T(1) - ratio) + T(N_Max))),
        N_Min);
    if (n < pop.pop_size()) {
      std::nth_element(
          pop.begin(), pop.begin() + n, pop.begin() + pop.pop_size(),
          [](auto &x, auto &y) { return x.fitness() < y.fitness(); });
      pop.resize(n);
    }
  }
};

} // namespace sevobench::de_module
