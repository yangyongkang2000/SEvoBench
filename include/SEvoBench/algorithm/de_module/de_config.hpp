#pragma once

#include "../evolutionary_algorithm.hpp"
#include "de_archive.hpp"
#include "de_constraint.hpp"
#include "de_crossover.hpp"
#include "de_mutation.hpp"
#include "de_parameter.hpp"
#include "de_population.hpp"
#include <memory>

namespace sevobench::de_module {
template <bool Use_Archive, std::floating_point T> struct de_config {
  std::unique_ptr<de_parameter<T>> parameter;
  std::unique_ptr<de_mutation<T>> mutation;
  std::unique_ptr<de_constraint<T>> constraint_handler;
  std::unique_ptr<de_crossover<T>> crossover;
  std::unique_ptr<de_population<T>> population_strategy;
  [[maybe_unused]] std::conditional_t<
      Use_Archive, std::unique_ptr<de_archive<T>>, void *> archive;

  void iterator(auto &&pop, auto &&trial, auto &&f, T lb, T ub,
                evolutionary_algorithm &alg) noexcept {
    parameter->prepare(pop);
    mutation->prepare(pop);
    crossover->prepare(pop);
    for (int i = 0; i < pop.pop_size(); i++) {
      do {
        if constexpr (Use_Archive) {
          mutation->mutate(trial[i], pop, archive->get(), parameter->get_f(i),
                           i);
        } else {
          mutation->mutate(trial[i], pop, parameter->get_f(i), i);
        }
        constraint_handler->repair(trial[i], lb, ub,
                                   pop[mutation->base_index(i)], pop[i]);
      } while (!constraint_handler->is_repair(
          trial[i], lb, ub, pop[mutation->base_index(i)], pop[i]));
      crossover->crossover(trial[i], pop[i], parameter->get_cr(i));
      trial[i].evaluate(f);
      alg.increment_fes();
    }
    parameter->update(pop, trial);
    for (int i = 0; i < pop.pop_size(); i++)
      if (trial[i].fitness() < pop[i].fitness()) {
        std::swap(trial[i], pop[i]);
        if constexpr (Use_Archive) {
          archive->add(trial[i]);
        }
      }
    population_strategy->resize(pop);
    alg.set_pop_size(pop.pop_size());
    if constexpr (Use_Archive) {
      archive->resize(pop);
    }
    alg.increment_iterator();
  }
};
} // namespace sevobench::de_module
