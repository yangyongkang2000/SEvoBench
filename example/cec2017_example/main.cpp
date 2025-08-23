#include "SEvoBench/sevobench.hpp"
#include <iostream>
int main() {
  using namespace sevobench::problem;

  // Create CEC2017 suite with F15 in 10D
  auto suite = suite_builder<cec2017>()
                   .dim<10>()
                   .type<float>()
                   .problem_index({15})
                   .instance_count(3)
                   .build();

  // Evaluate solution
  std::vector<float> x(10, 0.5f);
  for (auto &prob : suite) {
    auto fitness = (*prob)(x);
    std::cout << "F15 Instance " << prob->instance() << " fitness: " << fitness
              << "\n";
  }
}