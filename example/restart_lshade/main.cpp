#include <iostream>
#include "restart_lshade.h"  // Assuming the implementation is in this header

// Rosenbrock function implementation
template<std::floating_point T>
struct rosenbrock {
  T operator()(const T* x, int dim) const {
    T sum = 0;
    for (int i = 0; i < dim - 1; i++) {
      T t1 = x[i+1] - x[i]*x[i];
      T t2 = 1 - x[i];
      sum += 100 * t1*t1 + t2*t2;
    }
    return sum;
  }

  T lower_bound() const { return -100.0; }
  T upper_bound() const { return 100.0; }
};

int main() {
  const int dim = 30;            // Problem dimension
  const int max_fes = 10000*dim; // Maximum function evaluations
  const int pop_size = 18*dim;   // Population size (18*D as in the paper)
  const int B = 20;              // RL-SHADE parameter

  // Create RL-SHADE optimizer
  sevobench::rl_shade<double> optimizer(max_fes, dim, B);

  // Initialize population
  sevobench::population<double> pop(pop_size, dim,
                                    rosenbrock<double>{}.lower_bound(),
                                    rosenbrock<double>{}.upper_bound());

  // Run optimization
  optimizer.run(pop, rosenbrock<double>{});

  // Get results
  double best_fitness = optimizer.get_global_best_fitness();
  const auto& best_solution = optimizer.get_best_solution();

  std::cout << "Optimization completed!\n";
  std::cout << "Best fitness: " << best_fitness << "\n";
  std::cout << "Best solution: [";
  for (int i = 0; i < std::min(5, dim); i++) {  // Print first 5 dimensions
    std::cout << best_solution[i] << (i < dim-1 ? ", " : "");
  }
  if (dim > 5) std::cout << "...] (" << dim << " dimensions)";
  else std::cout << "]";

  return 0;
}