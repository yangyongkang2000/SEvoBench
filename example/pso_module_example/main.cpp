#include "SEvoBench/sevobench.hpp"
#include <iostream>
int main() {
  using namespace sevobench;
  using namespace pso_module;
  evolutionary_algorithm alg(10000, 50, 30); // max_fes=10^4
  // Configure SPSO
  auto builder =
      pso_algorithm_builder<double>()
          .update(std::make_unique<decrease_inertia_weight_update<double>>(
              alg, 0.4, 0.9, 1.494, 1.494))
          .topology(std::make_unique<gbest_topology<double>>())
          .constraint_handler(pso_constraint<double>{
              std::make_unique<spso_velocity_constraint<double>>(-2.0, 2.0),
              std::make_unique<positions_constraint<double>>()});

  auto pso = builder.build();

  // Initialize parameters
  population<double> pop(50, 30, -5.0, 5.0);
  pso_module::pso_velocity<double> vec(50, std::vector<double>(30));
  auto sphere_function = [](auto &x) {
    return std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  };
  // Run optimization
  pso.run(pop, vec, sphere_function, -5.0, 5.0, alg);

  // Output results
  std::cout << "Best fitness: " << pso.topology()->best_value() << std::endl;
}