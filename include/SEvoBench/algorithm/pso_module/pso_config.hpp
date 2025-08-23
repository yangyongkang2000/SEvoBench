#pragma once
#include "pso_constraint.hpp"
#include "pso_topology.hpp"
#include "pso_update.hpp"
namespace sevobench::pso_module {
template <std::floating_point T> struct pso_config {
  std::unique_ptr<pso_update<T>> update;
  std::unique_ptr<pso_topology<T>> topology;
  pso_constraint<T> constraint_handler;
  void iterator(auto &&pop, auto &&vec, auto &&f, T lb, T ub,
                evolutionary_algorithm &alg) noexcept {
    for (int i = 0; i < pop.pop_size(); i++) {
      update->update_velocity(pop, vec, *topology, i);
      constraint_handler.repair_velocity(vec[i]);
      update->update_position(pop, vec, *topology, i);
      constraint_handler.repair_position(pop[i], lb, ub);
      pop[i].evaluate(f);
      alg.increment_fes();
    }
    update->update(pop, *topology);
    topology->update(pop);
    alg.increment_iterator();
  }
};
} // namespace sevobench::pso_module
