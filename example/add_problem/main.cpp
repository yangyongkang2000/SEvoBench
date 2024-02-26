#include "SEvoBench/parallel_algorithm_benchmark.hpp"
#include "add_problem.hpp"
#include"fast_add_problem.hpp"
#include <iostream>
int main() {
  auto result =
      sevobench::evo_bench<sevobench::SLPSO, AddProblem, false, 40, 104>();
  for (auto _ : result.get_table_data())
    std::cout << _ << "\n";
  return 0;
}
