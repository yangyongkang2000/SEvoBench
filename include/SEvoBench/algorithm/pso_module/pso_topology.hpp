#pragma once

#include "../../common/population.hpp"
#include "../../common/tool.hpp"
#include "../evolutionary_algorithm.hpp"
#include "pso_velocity.hpp"
#include <algorithm>
#include <concepts>
#include <numeric>

namespace sevobench::pso_module {
template <std::floating_point T> class pso_topology {
public:
  virtual void prepare(const population<T> &) {};

  virtual solution<T> &local_best(int) = 0;

  virtual solution<T> &personal_best(int) = 0;

  virtual void update(const population<T> &) {};

  virtual solution<T> &best_solution() = 0;

  virtual T best_value() = 0;

  virtual int *neighbor(int, int *) = 0;

  virtual int *informant(int, int *) = 0;

  virtual ~pso_topology() {};
};
namespace detail {
template <std::floating_point T> class pbest_base : public pso_topology<T> {
protected:
  population<T> pbest;

public:
  solution<T> &personal_best(int i) override { return pbest[i]; }

  solution<T> &best_solution() override {
    return *std::min_element(pbest.begin(), pbest.end(), [](auto &x, auto &y) {
      return x.fitness() < y.fitness();
    });
  }

  T best_value() override {
    return std::min_element(
               pbest.begin(), pbest.end(),
               [](auto &x, auto &y) { return x.fitness() < y.fitness(); })
        ->fitness();
  }

  void update(const population<T> &pop) override {
    for (int i = 0; i < pop.pop_size(); i++) {
      if (pop[i].fitness() < pbest[i].fitness()) {
        pbest[i] = pop[i];
      }
    }
  }
};
} // namespace detail

template <std::floating_point T>
class gbest_topology final : public detail::pbest_base<T> {
  int min_index;

public:
  void prepare(const population<T> &pop) override {
    this->pbest = pop;
    min_index = static_cast<int>(
        std::min_element(
            this->pbest.begin(), this->pbest.end(),
            [](auto &l, auto &r) { return l.fitness() < r.fitness(); }) -
        this->pbest.begin());
  }

  solution<T> &local_best(int) override { return this->pbest[min_index]; }

  solution<T> &best_solution() override { return this->pbest[min_index]; }

  T best_value() override { return this->pbest[min_index].fitness(); }

  void update(const population<T> &pop) override {
    for (int i = 0; i < pop.pop_size(); i++) {
      if (pop[i].fitness() < this->pbest[i].fitness()) {
        this->pbest[i] = pop[i];
      }
    }
    min_index = static_cast<int>(
        std::min_element(
            this->pbest.begin(), this->pbest.end(),
            [](auto &l, auto &r) { return l.fitness() < r.fitness(); }) -
        this->pbest.begin());
  }

  int *neighbor(int, int *o) override {
    std::iota(o, o + this->pbest.pop_size(), int(0));
    return o + this->pbest.pop_size();
  }

  int *informant(int i, int *o) override {
    *o++ = i;
    if (i != min_index)
      *o++ = min_index;
    return o;
  }
};
namespace detail {
template <int _K, std::floating_point T>
class kneighbor_topology : public pbest_base<T> {
protected:
  std::vector<std::array<int, _K>> topology_index;
  static constexpr int K = _K;

public:
  solution<T> &local_best(int i) override {
    auto index = *std::min_element(
        topology_index[i].begin(), topology_index[i].end(),
        [=, this](auto i1, auto i2) {
          return this->pbest[i1].fitness() < this->pbest[i2].fitness();
        });
    return this->pbest[index];
  }

  int *neighbor(int i, int *o) override {
    std::copy_n(topology_index[i].data(), _K, o);
    std::sort(o, o + _K);
    return std::unique(o, o + _K);
  }

  int *informant(int i, int *o) override {
    *o++ = i;
    auto index = *std::min_element(
        topology_index[i].begin(), topology_index[i].end(),
        [=, this](auto i1, auto i2) {
          return this->pbest[i1].fitness() < this->pbest[i2].fitness();
        });
    if (i != index) {
      *o++ = index;
    }
    return o;
  }
};
} // namespace detail
template <std::floating_point T>
class lbest_topology final : public detail::kneighbor_topology<3, T> {
public:
  void prepare(const population<T> &pop) override {
    this->pbest = pop;
    this->topology_index.resize(pop.pop_size());
    for (int i = 0; i < pop.pop_size(); i++) {
      for (int j = 0; j < this->K; j++)
        this->topology_index[i][j] = (i + j - 1) % pop.pop_size();
      this->topology_index[0][0] = (pop.pop_size() - 1);
    }
  }
};

template <std::floating_point T, int K = 4, typename R = sevobench::tool::rng>
  requires tool::random_generator_concept<R, T>

class random_topology final : public detail::kneighbor_topology<K, T> {
  const int m = 0;
  int n = 0;

public:
  R RNG;

  random_topology() = default;

  random_topology(int _m) : m(_m) {}

  void prepare(const population<T> &pop) override {
    this->pbest = pop;
    this->topology_index.resize(pop.pop_size());
    for (int i = 0; i < pop.pop_size(); i++) {
      this->topology_index[i][0] = i;
      for (int j = 1; j < this->K; j++)
        this->topology_index[i][j] = RNG.rand_int(pop.pop_size());
    }
  }

  void update(const population<T> &pop) override {
    auto old_value = this->best_value();
    for (int i = 0; i < pop.pop_size(); i++) {
      if (pop[i].fitness() < this->pbest[i].fitness()) {
        this->pbest[i] = pop[i];
      }
    }
    auto new_value = this->best_value();
    if (new_value < old_value) {
      n = 0;
    } else {
      if (++n > m) {
        for (int i = 0; i < pop.pop_size(); i++) {
          this->topology_index[i][0] = i;
          for (int j = 1; j < this->K; j++)
            this->topology_index[i][j] = RNG.rand_int(pop.pop_size());
        }
        n = 0;
      }
    }
  }
};

template <std::floating_point T, typename R = sevobench::tool::rng>
  requires tool::random_generator_concept<R, T>

class von_neumann_topology final : public detail::kneighbor_topology<5, T> {
public:
  void prepare(const population<T> &pop) override {
    this->pbest = pop;
    this->topology_index.resize(pop.pop_size());
    int rows = static_cast<int>(std::sqrt(pop.pop_size()));
    while (pop.pop_size() % rows != 0) {
      rows--;
    }
    int const columns = pop.pop_size() / rows;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        this->topology_index[i * columns + j][4] = i * columns + j;
        if (i != rows - 1)
          this->topology_index[i * columns + j][0] = (i + 1) * columns + j;
        else
          this->topology_index[i * columns + j][0] = j;

        if (i != 0)
          this->topology_index[i * columns + j][1] = (i - 1) * columns + j;
        else
          this->topology_index[i * columns + j][1] = (rows - 1) * columns + j;

        if (j != columns - 1)
          this->topology_index[i * columns + j][2] = i * columns + (j + 1);
        else
          this->topology_index[i * columns + j][2] = i * columns;

        if (j != 0)
          this->topology_index[i * columns + j][3] = i * columns + (j - 1);
        else
          this->topology_index[i * columns + j][3] =
              i * columns + (columns - 1);
      }
    }
  }
};

template <bool b, std::floating_point T, typename R = sevobench::tool::rng>
  requires tool::random_generator_concept<R, T>

class time_varying_topology : public detail::pbest_base<T> {
  std::vector<std::vector<int>> topology_index;
  int current_iterator;
  int per_iterator;
  int current_size;

public:
  R RNG;

  time_varying_topology(int _per_iterator) : per_iterator(_per_iterator) {}

  void prepare(const population<T> &pop) override {
    this->pbest = pop;
    topology_index.resize(pop.pop_size());
    for (int i = 0; i < pop.pop_size(); i++) {
      topology_index[i].resize(pop.pop_size());
      std::iota(topology_index[i].begin(), topology_index[i].end(), 0);
      std::partition(topology_index[i].begin(), topology_index[i].end(),
                     [i, len = pop.pop_size()](auto x) {
                       if (x == i || x == (i + 1) % len ||
                           x == ((len + i - 1) % len))
                         return true;
                       return false;
                     });
    }
    current_size = b ? pop.pop_size() : 3;
    current_iterator = 0;
  };

  solution<T> &local_best(int i) override {
    auto index = *std::min_element(
        topology_index[i].begin(), topology_index[i].begin() + current_size,
        [=, this](auto i1, auto i2) {
          return this->pbest[i1].fitness() < this->pbest[i2].fitness();
        });
    return this->pbest[index];
  }

  void update(const population<T> &pop) override {
    detail::pbest_base<T>::update(pop);
    if (++current_iterator % per_iterator == 0) {
      if constexpr (b) {
        if (current_size > 3) {
          for (auto &_ : topology_index) {
            auto idx = RNG.rand_int(current_size - 3);
            std::swap(_[idx + 3], _.back());
            _.pop_back();
          }
          --current_size;
        }
      } else {
        if (current_size < pop.pop_size()) {
          for (auto &_ : topology_index) {
            auto idx = RNG.rand_int(pop.pop_size() - current_size);
            std::swap(_[idx + current_size], _[current_size]);
          }
          ++current_size;
        }
      }
    }
  }

  int *neighbor(int i, int *o) override {
    return std::copy_n(topology_index[i].data(), current_size, o);
  }

  int *informant(int i, int *o) override {
    *o++ = i;
    auto index = *std::min_element(
        topology_index[i].begin(), topology_index[i].begin() + current_size,
        [=, this](auto i1, auto i2) {
          return this->pbest[i1].fitness() < this->pbest[i2].fitness();
        });
    if (i != index) {
      *o++ = index;
    }
    return o;
  }
};

template <std::floating_point T, typename R = sevobench::tool::rng>
using decreasing_topology = time_varying_topology<true, T, R>;
template <std::floating_point T, typename R = sevobench::tool::rng>
using increasing_topology = time_varying_topology<false, T, R>;

template <std::floating_point T, int K = 3, typename R = sevobench::tool::rng>
  requires tool::random_generator_concept<R, T>

class dms_topology final : public detail::kneighbor_topology<K, T> {
  std::vector<int> indices;
  const int regroup_interval;
  int current_iter;

public:
  R RNG;

  dms_topology(int interval) : regroup_interval(interval) {}

  void prepare(const population<T> &pop) override {
    this->pbest = pop;
    current_iter = 0;
    this->topology_index.resize(pop.pop_size());
    indices.resize(pop.pop_size());
    std::iota(indices.begin(), indices.end(), 0);
    regroup_subswarms();
  }

  void regroup_subswarms() {
    tool::kunth_shuffle(indices.begin(), indices.end(), RNG());
    const int size = static_cast<int>(indices.size());
    for (int i = 0; i < size; i += K) {
      const int K1 = i + K < size ? K : size - i;
      for (int j = 0; j < K1; j++) {
        auto &neighbors = this->topology_index[indices[i + j]];
        neighbors[0] = indices[i + j];
        neighbors[1] = indices[i + (j + 1) % K1];
        neighbors[2] = indices[i + (j + K1 - 1) % K1];
      }
    }
  }

  void update(const population<T> &pop) override {
    detail::kneighbor_topology<K, T>::update(pop);
    if (++current_iter % regroup_interval == 0)
      regroup_subswarms();
  }
};

} // namespace sevobench::pso_module
