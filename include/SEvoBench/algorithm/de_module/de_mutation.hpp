
#pragma once

#include "../../common/population.hpp"
#include "../evolutionary_algorithm.hpp"

namespace sevobench::de_module {

template <std::floating_point T> class de_mutation {
public:
  virtual void prepare(population<T> &) {};

  virtual void mutate(solution<T> &, const population<T> &, T, int) = 0;

  virtual void mutate(solution<T> &, const population<T> &,
                      std::span<const solution<T>>, T, int) {};

  virtual ~de_mutation() = default;

  virtual int base_index(int i) { return i; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class rand1_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r1, r2, r3] = RNG.template pick_random<3>(pop.pop_size(), i);
    index = r1;
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[r1][j] + f * (pop[r2][j] - pop[r3][j]);
    }
  }

  int base_index(int) override { return index; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class rand2_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r1, r2, r3, r4, r5] = RNG.template pick_random<5>(pop.pop_size(), i);
    index = r1;
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] =
          pop[r1][j] + f * (pop[r2][j] - pop[r3][j] + pop[r4][j] - pop[r5][j]);
    }
  }

  int base_index(int) override { return index; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class best1_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void prepare(population<T> &pop) override {
    index =
        static_cast<int>(std::min_element(pop.begin(), pop.end(),
                                          [](auto &x, auto &y) {
                                            return x.fitness() < y.fitness();
                                          }) -
                         pop.begin());
  }

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r1, r2] = RNG.template pick_random<2>(pop.pop_size(), i);
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[index][j] + f * (pop[r1][j] - pop[r2][j]);
    }
  }

  int base_index(int) override { return index; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class best2_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void prepare(population<T> &pop) override {
    index =
        static_cast<int>(std::min_element(pop.begin(), pop.end(),
                                          [](auto &x, auto &y) {
                                            return x.fitness() < y.fitness();
                                          }) -
                         pop.begin());
  }

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r1, r2, r3, r4] = RNG.template pick_random<4>(pop.pop_size(), i);
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[index][j] +
                 f * (pop[r1][j] - pop[r2][j] + pop[r3][j] - pop[r4][j]);
    }
  }

  int base_index(int) override { return index; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class ttpb1_mutation final : public de_mutation<T> {
  const T p;
  int P;

public:
  R RNG;

  ttpb1_mutation(T _p = T(0.11)) : p(_p) {}

  void prepare(population<T> &pop) override {
    P = static_cast<int>(std::max(pop.pop_size() * p, T(2)));
    std::nth_element(
        pop.begin(), pop.begin() + P, pop.begin() + pop.pop_size(),
        [](auto &x, auto &y) { return x.fitness() < y.fitness(); });
  }

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto r0 = RNG.rand_int(P);
    auto [r1, r2] = RNG.template pick_random<2>(pop.pop_size(), i);
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] =
          pop[i][j] + f * (pop[r0][j] - pop[i][j] + pop[r1][j] - pop[r2][j]);
    }
  }

  void mutate(solution<T> &trial, const population<T> &pop,
              std::span<const solution<T>> archives, T f, int i) override {
    auto r0 = RNG.rand_int(P);
    int r1, r2;
    do {
      r1 = RNG.rand_int(pop.pop_size());
    } while (r1 == i);
    do {
      r2 = RNG.rand_int(pop.pop_size() + static_cast<int>(archives.size()));
    } while (r2 == i || r2 == r1);
    const auto &tmp =
        r2 >= pop.pop_size() ? archives[r2 - pop.pop_size()] : pop[r2];
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[i][j] + f * (pop[r0][j] - pop[i][j] + pop[r1][j] - tmp[j]);
    }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class ttb1_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void prepare(population<T> &pop) override {
    index = int(std::min_element(pop.begin(), pop.end(),
                                 [](auto &x, auto &y) {
                                   return x.fitness() < y.fitness();
                                 }) -
                pop.begin());
  }

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r1, r2] = RNG.template pick_random<2>(pop.pop_size(), i);
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] =
          pop[i][j] + f * (pop[index][j] - pop[i][j] + pop[r1][j] - pop[r2][j]);
    }
  }

  void mutate(solution<T> &trial, const population<T> &pop,
              std::span<const solution<T>> archives, T f, int i) override {
    int r1, r2;
    do {
      r1 = RNG.rand_int(pop.pop_size());
    } while (r1 == i);
    do {
      r2 = RNG.rand_int(pop.pop_size() + static_cast<int>(archives.size()));
    } while (r2 == i || r2 == r1);
    const auto &tmp =
        r2 >= pop.pop_size() ? archives[r2 - pop.pop_size()] : pop[r2];
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] =
          pop[i][j] + f * (pop[index][j] - pop[i][j] + pop[r1][j] - tmp[j]);
    }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class ttb2_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void prepare(population<T> &pop) override {
    index = int(std::min_element(pop.begin(), pop.end(),
                                 [](auto &x, auto &y) {
                                   return x.fitness() < y.fitness();
                                 }) -
                pop.begin());
  }

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r1, r2, r3, r4] = RNG.template pick_random<4>(pop.pop_size(), i);
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[i][j] + f * (pop[index][j] - pop[i][j] + pop[r1][j] -
                                  pop[r2][j] + pop[r3][j] - pop[r4][j]);
    }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class two_opt1_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r0, r1, r2] = RNG.template pick_random<3>(pop.pop_size(), i);
    if (pop[r0].fitness() > pop[r1].fitness())
      std::swap(r0, r1);
    index = r0;
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[r0][j] + f * (pop[r1][j] - pop[r2][j]);
    }
  }

  int base_index(int) override { return index; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class two_opt2_mutation final : public de_mutation<T> {
  int index;

public:
  R RNG;

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto [r0, r1, r2, r3, r4] = RNG.template pick_random<5>(pop.pop_size(), i);
    if (pop[r0].fitness() > pop[r1].fitness())
      std::swap(r0, r1);
    index = r0;
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] =
          pop[r0][j] + f * (pop[r1][j] - pop[r2][j] + pop[r3][j] - pop[r4][j]);
    }
  }

  int base_index(int) override { return index; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class ttpb1_weight_mutation final : public de_mutation<T> {
  const evolutionary_algorithm &alg;
  const T p_min;
  const T p_max;
  int P;

public:
  R RNG;

  ttpb1_weight_mutation(const evolutionary_algorithm &_alg, T _p_min, T _p_max)
      : alg(_alg), p_min(_p_min), p_max(_p_max) {}

  void prepare(population<T> &pop) override {
    auto p = (p_max - p_min) / alg.max_fes() * alg.current_fes() + p_min;
    P = static_cast<int>(std::max(pop.pop_size() * p, T(2)));
    std::nth_element(
        pop.begin(), pop.begin() + P, pop.begin() + pop.pop_size(),
        [](auto &x, auto &y) { return x.fitness() < y.fitness(); });
  }

  void mutate(solution<T> &trial, const population<T> &pop, T f,
              int i) override {
    auto weight =
        alg.current_fes() < alg.max_fes() / 5
            ? T(0.7)
            : (alg.current_fes() < (2 * alg.max_fes() / 5) ? T(0.8) : T(1.2));
    auto r0 = RNG.rand_int(P);
    auto [r1, r2] = RNG.template pick_random<2>(pop.pop_size(), i);
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[i][j] + f * (weight * (pop[r0][j] - pop[i][j]) +
                                  pop[r1][j] - pop[r2][j]);
    }
  }

  void mutate(solution<T> &trial, const population<T> &pop,
              std::span<const solution<T>> archives, T f, int i) override {
    auto weight =
        alg.current_fes() < alg.max_fes() / 5
            ? T(0.7)
            : (alg.current_fes() < (2 * alg.max_fes() / 5) ? T(0.8) : T(1.2));
    auto r0 = RNG.rand_int(P);
    int r1, r2;
    do {
      r1 = RNG.rand_int(pop.pop_size());
    } while (r1 == i);
    do {
      r2 = RNG.rand_int(pop.pop_size() + static_cast<int>(archives.size()));
    } while (r2 == i || r2 == r1);
    const auto &tmp =
        r2 >= pop.pop_size() ? archives[r2 - pop.pop_size()] : pop[r2];
    for (int j = 0; j < pop.dim(); j++) {
      trial[j] = pop[i][j] +
                 f * (weight * (pop[r0][j] - pop[i][j]) + pop[r1][j] - tmp[j]);
    }
  }
};

} // namespace sevobench::de_module
