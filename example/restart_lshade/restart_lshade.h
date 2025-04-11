#pragma once
#include "SEvoBench/sevobench.hpp"
#include <algorithm>
#include <memory>
#include <numeric>

namespace sevobench {

// RL-SHADE algorithm controller that manages restart conditions and tracks
// global best fitness
template <std::floating_point T = float>
class rl_shade_algorithm : public evolutionary_algorithm {
protected:
  T global_best_fitness; // Stores the best fitness value found during
                         // optimization
  const int max_evals_per_iter_base; // Base value for max evaluations per
                                     // iteration (MaxFEs/B)

public:
  const int improvement_threshold; // Threshold for restart (500 * dimension)

  // Constructor initializes algorithm parameters
  rl_shade_algorithm(int max_fes, int dim, int B = 10)
      : evolutionary_algorithm(max_fes),
        global_best_fitness(
            std::numeric_limits<T>::max()), // Initialize to worst possible
                                            // value
        max_evals_per_iter_base(max_fes /
                                B), // Calculate base evaluations per iteration
        improvement_threshold(500 * dim) { // Set restart threshold
    set_dim(dim);                          // Initialize problem dimension
  }

  // Returns the maximum evaluations for current iteration (dynamically
  // adjusted)
  int get_max_evals_per_iter() const {
    return std::min(max_evals_per_iter_base,
                    max_fes() -
                        current_fes()); // Ensure we don't exceed total budget
  }

  // Updates the global best fitness if a better solution is found
  void update_best_fitness(T current_fitness) {
    if (current_fitness < global_best_fitness) {
      global_best_fitness = current_fitness;
    }
  }

  // Returns the current global best fitness value
  T get_global_best() const { return global_best_fitness; }
};

// Main RL-SHADE implementation class
template <std::floating_point T = float> class rl_shade {
  rl_shade_algorithm<T> alg; // Algorithm state and restart controller
  de_module::de_algorithm<true, T> shade; // SHADE optimization core
  solution<T> best_solution;              // Stores the best solution found

public:
  // Constructor sets up RL-SHADE components with specified parameters
  rl_shade(int max_fes, int dim, int B = 10)
      : alg(max_fes, dim, B), // Initialize algorithm controller
        shade(
            de_module::de_algorithm_builder<T>()
                // Mutation strategy: current-to-pbest/1
                .mutation(std::make_unique<de_module::ttb1_mutation<T>>())
                // Boundary constraint handling
                .constraint_handler(
                    std::make_unique<de_module::rand_base_repair<T>>())
                // Crossover operator
                .crossover(std::make_unique<de_module::binomial_crossover<T>>())
                // Parameter adaptation
                .parameter(std::make_unique<de_module::shade_parameter<T>>())
                // Population size reduction strategy
                .population_strategy(
                    std::make_unique<de_module::linear_reduction<T>>(alg, 4,
                                                                     18 * dim))
                // External archive implementation
                .archive(std::make_unique<de_module::fifo_archive<T>>())
                .build()),
        best_solution(dim) // Initialize the best solution storage
  {}                       // Build the SHADE core

  // Main optimization routine with explicit bounds
  void run(population<T> &pop, auto &&f, T lb, T ub) {
    population<T> trial(pop.pop_size(),
                        pop.dim()); // Initialize trial population

    // Main optimization loop - runs until evaluation budget exhausted
    while (alg.current_fes() < alg.max_fes()) {
      const int current_max_evals =
          alg.get_max_evals_per_iter(); // Get current iteration budget

      // Evaluate initial population
      for (auto &_ : pop)
        _.evaluate(f);
      alg.add_fes(pop.pop_size()); // Update evaluation count

      // Track the best solution at iteration start
      T iteration_start_best =
          std::min_element(pop.begin(), pop.end(),
                           [](const auto &x, const auto &y) {
                             return x.fitness() < y.fitness();
                           })
              ->fitness();
      int start_fes = alg.current_fes(); // Track starting evaluation count
      const int least_stop_fes =
          start_fes + current_max_evals; // Minimum stopping point
      shade.archive()->prepare(pop);     // Prepare archive

      // Inner optimization loop
      do {
        // Generate and evaluate trial solutions
        shade.iterator(pop, trial, f, lb, ub, alg);
      } while ([&, least_stop_fes](const auto &current_pop, const auto &a) {
        // Find current best solution
        auto current_best =
            std::min_element(current_pop.begin(), current_pop.end(),
                             [](const auto &x, const auto &y) {
                               return x.fitness() < y.fitness();
                             });

        // Update global best if improved
        if (current_best->fitness() < alg.get_global_best()) {
          best_solution = *current_best;
          alg.update_best_fitness(current_best->fitness());
        }

        // Check restart conditions:
        // 1. No improvement in this iteration
        bool no_improvement = (current_best->fitness() >= iteration_start_best);
        iteration_start_best = current_best->fitness();
        int evals_used = a.current_fes() - start_fes;

        // Reset counter if improvement found
        if (!no_improvement) {
          start_fes = a.current_fes();
        }

        // Continue optimization if:
        // 1. Haven't reached minimum evaluations OR
        // 2. Haven't exceeded max FEs AND
        // 3. Haven't reached improvement threshold
        if (a.current_fes() < least_stop_fes) {
          return true;
        }
        if (a.current_fes() > a.max_fes()) {
          return false;
        }
        return evals_used < alg.improvement_threshold;
      }(pop, alg));

      // Restart procedure:
      // 1. Reinitialize population randomly
      // 2. Reset parameter adaptation memory
      pop.randomize(lb, ub);
      static_cast<de_module::shade_parameter<T> *>(shade.parameter())->reset();
    }
  }

  // Bound-constrained optimization using function's native bounds
  void run(population<T> &pop, auto &&f) {
    run(pop, f, f.lower_bound(), f.upper_bound());
  }

  // Returns best fitness value found
  T get_global_best_fitness() const { return alg.get_global_best(); }

  // Returns reference to best solution found
  const solution<T> &get_best_solution() const { return best_solution; }
};

} // namespace sevobench