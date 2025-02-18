#pragma once
#include "pso_config.hpp"
namespace sevobench::pso_module {
template <std::floating_point T> class pso_algorithm {
  pso_config<T> config;

public:
  using value_type = T;
  pso_algorithm(pso_config<T> &&_config) : config(std::move(_config)) {}
  template <typename F = evolutionary_algorithm_condition>
  void run(std::convertible_to<population<T>> auto &&pop,
           std::convertible_to<pso_velocity<T>> auto &&vec, auto &&f, T lb,
           T ub, evolutionary_algorithm &alg, F condition = {}) noexcept {
    alg.set_pop_size(pop.pop_size());
    alg.set_dim(pop.dim());
    for (auto &_ : pop)
      _.evaluate(f);
    alg.add_fes(pop.pop_size());
    config.topology->prepare(pop);
    do {
      config.iterator(pop, vec, f, lb, ub, alg);
    } while (condition(pop, alg));
  }
  template <typename F = evolutionary_algorithm_condition>
  void run(std::convertible_to<population<T>> auto &&pop,
           std::convertible_to<pso_velocity<T>> auto &&vec, auto &&f,
           evolutionary_algorithm &alg, F condition = {}) noexcept {
    run(pop, vec, f, f.lower_bound(), f.upper_bound(), alg, condition);
  }
  auto update() noexcept { return config.update.get(); }
  auto replace_update(std::unique_ptr<pso_update<T>> &&ptr) noexcept {
    config.update = std::move(ptr);
  }
  auto topology() noexcept { return config.topology.get(); }
  auto replace_topology(std::unique_ptr<pso_topology<T>> &&ptr) noexcept {
    config.topology = std::move(ptr);
  }
  auto &constraint_handler() noexcept { return config.constraint_handler; }
  auto replace_constraint_handler(pso_constraint<T> &&pc) noexcept {
    config.constraint_handler = std::move(pc);
  }
};

template <std::floating_point T = float, bool B1 = false, bool B2 = false>
class [[nodiscard]] pso_algorithm_builder {
  pso_config<T> config_;

public:
  pso_algorithm_builder() = default;
  pso_algorithm_builder(pso_config<T> &&_config)
      : config_(std::move(_config)) {}
  template <std::floating_point T1> auto type() noexcept {
    return pso_algorithm_builder<T1, B1, B2>(std::move(config_));
  }
  auto update(std::unique_ptr<pso_update<T>> &&ptr) noexcept {
    config_.update = std::move(ptr);
    return pso_algorithm_builder<T, true, B2>(std::move(config_));
  }
  auto topology(std::unique_ptr<pso_topology<T>> &&ptr) noexcept {
    config_.topology = std::move(ptr);
    return pso_algorithm_builder<T, B1, true>(std::move(config_));
  }
  auto constraint_handler(pso_constraint<T> &&pc) noexcept {
    config_.constraint_handler = std::move(pc);
    return pso_algorithm_builder<T, B1, B2>(std::move(config_));
  }
  auto build() noexcept {
    static_assert(B1, "PSO UPDATE IS MISSING!");
    static_assert(B2, "PSO TOPOLOGY IS MISSING!");
    return pso_algorithm<T>(std::move(config_));
  }
};

} // namespace sevobench::pso_module
