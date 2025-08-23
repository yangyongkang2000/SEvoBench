
#pragma once

#include "../../common/population.hpp"
#include "../../common/tool.hpp"

namespace sevobench::de_module {
template <std::floating_point T> class de_crossover {
public:
  virtual void crossover(std::span<T>, std::span<const T>, T) = 0;

  virtual void prepare(const population<T> &) {};

  virtual ~de_crossover() {};
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class binomial_crossover final : public de_crossover<T> {
public:
  R RNG;

  void crossover(std::span<T> donor, std::span<const T> target, T cr) override {
    auto dim = static_cast<int>(target.size());
    auto j = RNG.rand_int(dim);
    for (int i = 0; i < dim; i++)
      donor[i] =
          (RNG.template rand_float<T>() < cr || i == j) ? donor[i] : target[i];
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class exponential_crossover final : public de_crossover<T> {
public:
  R RNG;

  void crossover(std::span<T> donor, std::span<const T> target, T cr) override {
    auto dim = static_cast<int>(target.size());
    int start = RNG.rand_int(dim);
    int L = 0;
    while (L < dim && RNG.template rand_float<T>() < cr)
      L++;
    int end = (start + L - 1) % dim;
    if (end < start) {
      for (int i = end + 1; i < start; i++) {
        donor[i] = target[i];
      }
    } else {
      for (int i = 0; i < start; i++) {
        donor[i] = target[i];
      }
      for (int i = end + 1; i < dim; i++) {
        donor[i] = target[i];
      }
    }
  }
};

    template <std::floating_point T, typename R = tool::rng>
    requires tool::random_generator_concept<R, T>

    class cauchy_pertubation final : public de_crossover<T>
    {
        const T pr=T(0.2);
        const T lb;
        const T ub;
    public:
        R RNG;
        R RNG1;
        R RNG2;
        cauchy_pertubation(T _lb,T _ub,T _pr=T(0.2)):lb(_lb),ub(_ub),pr(_pr) {}
        void crossover(std::span<T> donor, std::span<const T> target, T cr) override {
            auto dim = static_cast<int>(target.size());
            auto j = RNG.rand_int(dim);
            for (int i = 0; i < dim; i++)
                donor[i] =
                        (RNG.template rand_float<T>() < cr || i == j) ? donor[i] :(RNG1.template rand_float<T>()<pr?std::clamp(RNG2.cauchy(target[i],T(0.1)),lb,ub) :target[i]);
        }
    };
} // namespace sevobench::de_module
