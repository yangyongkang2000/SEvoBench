#include "SEvoBench/sevobench.hpp"
#include<iostream>
int main() {
  using namespace sevobench;
  using namespace de_module;
  // Define Sphere function
  auto sphere = [](auto &x) {
    return std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  };
  evolutionary_algorithm alg(10000, 50, 30); // max_fes=10^4
  // Configure JADE
  de_config<true, double> config{
      std::make_unique<jade_parameter<double>>(0.1, 0.5),
      std::make_unique<ttpb1_mutation<double>>(0.1),
      std::make_unique<projection_repair<double>>(),
      std::make_unique<binomial_crossover<double>>(),
      std::make_unique<linear_reduction<double>>(alg, 50, 100),
      std::make_unique<fifo_archive<double>>(2.0)
  };

  // Initialize algorithm
  de_algorithm<true, double> de(std::move(config));

  // Run optimization
  population<double> pop(50, 30, -5.0, 5.0);
  de.run(pop, sphere, -5.0, 5.0, alg);

  // Output results
  auto best = *std::min_element(pop.begin(), pop.end());
  std::cout << "Best fitness: " << best.fitness() << std::endl;
}