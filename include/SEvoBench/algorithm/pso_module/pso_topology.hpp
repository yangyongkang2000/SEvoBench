#pragma once
#include "../../common/population.hpp"
#include "../../common/tool.hpp"
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

} // namespace sevobench::pso_module
