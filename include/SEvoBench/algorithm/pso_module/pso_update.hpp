#pragma once

#include "../../algorithm/evolutionary_algorithm.hpp"
#include "../../common/population.hpp"
#include "../../common/tool.hpp"
#include "pso_topology.hpp"
#include "pso_velocity.hpp"

namespace sevobench::pso_module {
template <std::floating_point T> class pso_update {
public:
  virtual void update_velocity(const population<T> &, pso_velocity<T> &,
                               pso_topology<T> &, int) {};

  virtual void update_position(population<T> &pop, const pso_velocity<T> &vec,
                               pso_topology<T> &, int i) {
    int len = pop.dim();
    for (int j = 0; j < len; j++)
      pop[i][j] += vec[i][j];
  }

  virtual void update(const population<T> &, pso_topology<T> &) {}

  virtual ~pso_update() = default;
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class decrease_inertia_weight_update final : public pso_update<T> {
  const evolutionary_algorithm &alg;
  const T w_min = T(0.4);
  const T w_max = T(0.9);
  const T c1 = T(2);
  const T c2 = T(2);

public:
  R RNG1;
  R RNG2;

  decrease_inertia_weight_update(const evolutionary_algorithm &_alg,
                                 T _w_min = T(0.4), T _w_max = T(0.9),
                                 T _c1 = T(2), T _c2 = T(2))
      : alg(_alg), w_min(_w_min), w_max(_w_max), c1(_c1), c2(_c2) {}

  void update_velocity(const population<T> &pop, pso_velocity<T> &vec,
                       pso_topology<T> &top, int i) override {
    int dim = pop.dim();
    auto &pbest = top.personal_best(i);
    auto &gbest = top.local_best(i);
    auto w =
        w_max + (w_min - w_max) * alg.current_iterator() / alg.max_iterator();
    for (int j = 0; j < dim; j++) {
      vec[i][j] =
          w * vec[i][j] +
          c1 * (RNG1.template rand_float<T>()) * (pbest[j] - pop[i][j]) +
          c2 * (RNG2.template rand_float<T>()) * (gbest[j] - pop[i][j]);
    }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class inertia_weight_update final : public pso_update<T> {
  const T w = T(1) / T(2 * std::numbers::ln2_v<T>);
  const T c1 = T(0.5) + std::numbers::ln2_v<T>;
  const T c2 = T(0.5) + std::numbers::ln2_v<T>;

public:
  R RNG1;
  R RNG2;

  inertia_weight_update() = default;

  inertia_weight_update(T _w, T _c1, T _c2) : w(_w), c1(_c1), c2(_c2) {}

  void update_velocity(const population<T> &pop, pso_velocity<T> &vec,
                       pso_topology<T> &top, int i) override {
    int dim = pop.dim();
    auto &pbest = top.personal_best(i);
    auto &gbest = top.local_best(i);
    for (int j = 0; j < dim; j++) {
      vec[i][j] =
          w * vec[i][j] +
          c1 * (RNG1.template rand_float<T>()) * (pbest[j] - pop[i][j]) +
          c2 * (RNG2.template rand_float<T>()) * (gbest[j] - pop[i][j]);
    }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class bare_bone_update final : public pso_update<T> {
public:
  R RNG;

  void update_position(population<T> &pop, const pso_velocity<T> &,
                       pso_topology<T> &top, int i) override {
    int dim = pop.dim();
    auto &pbest = top.personal_best(i);
    auto &gbest = top.local_best(i);
    for (int j = 0; j < dim; j++) {
      pop[i][j] = RNG.normal(T(0.5) * (pbest[j] + gbest[j]),
                             std::abs(pbest[j] - gbest[j]));
    }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class fips_update : public pso_update<T> {
  const T phi = T(4.1);
  const T chi = T(0.729844);
  std::vector<int> index;

public:
  R RNG;

  fips_update() = default;

  fips_update(T _phi)
      : phi(_phi), chi(T(2) / (_phi - 2 + std::sqrt(_phi * (_phi - 4)))) {}

  void update_velocity(const population<T> &pop, pso_velocity<T> &vec,
                       pso_topology<T> &top, int i) override {
    int dim = pop.dim();
    if (static_cast<int>(index.size()) < pop.pop_size()) {
      index.resize(pop.pop_size());
    }
    std::span<int> s(index.data(), top.neighbor(i, index.data()));
    auto ratio = T(1) / T(s.size());
    for (auto v : s) {
      auto &pbest = top.personal_best(v);
      for (int j = 0; j < dim; j++)
        vec[i][j] += phi * (RNG.template rand_float<T>()) *
                     (pbest[j] - pop[i][j]) * ratio;
    }
    for (int j = 0; j < dim; j++)
      vec[i][j] *= chi;
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class spherical_update : public pso_update<T> {
  const T w = T(1) / T(2 * std::numbers::ln2_v<T>);
  const T c1 = T(0.5) + std::numbers::ln2_v<T>;
  const T c2 = T(0.5) + std::numbers::ln2_v<T>;
  std::vector<T> gpos;
  std::vector<T> xpos;

public:
  R RNG1;
  R RNG2;
  R RNG3;
  R RNG4;

  spherical_update() = default;

  spherical_update(T _w, T _c1, T _c2) : w(_w), c1(_c1), c2(_c2) {}

  void update_velocity(const population<T> &pop, pso_velocity<T> &vec,
                       pso_topology<T> &top, int i) override {
    auto hypersphere_sample = [](T *first, T r, T *out, int Dim) {
      auto sum = std::sqrt(std::inner_product(out, out + Dim, out, T(0)));
      for (int i = 0; i < Dim; i++)
        out[i] = first[i] + r * out[i] / sum;
      return out + Dim;
    };
    if (static_cast<int>(gpos.size()) < pop.dim()) {
      gpos.resize(pop.dim());
      xpos.resize(pop.dim());
    }
    T rsquare(0);
    auto &pbest = top.personal_best(i);
    auto &lbest = top.local_best(i);
    for (int j = 0; j < pop.dim(); j++) {
      auto t = pop[i][j];
      T tmp = (c1 * (RNG1.template rand_float<T>()) * (pbest[j] - t) +
               c2 * (RNG2.template rand_float<T>()) * (lbest[j] - t)) /
              T(3);
      rsquare += tmp * tmp;
      gpos[j] = tmp + t;
    }
    for (int j = 0; j < pop.dim(); j++)
      xpos[j] = RNG3.rand_float(T(-1), T(1));
    hypersphere_sample(gpos.data(),
                       std::sqrt(rsquare) * (RNG4.template rand_float<T>()),
                       xpos.data(), pop.dim());
    for (int j = 0; j < pop.dim(); j++)
      vec[i][j] = w * vec[i][j] + xpos[j] - pop[i][j];
  }
};

} // namespace sevobench::pso_module
