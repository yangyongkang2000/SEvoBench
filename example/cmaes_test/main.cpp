#include <iostream>

#include "SEvoBench/parallel_algorithm_benchmark.hpp"
#include "cmaes_test.hpp"
int main() {
  constexpr int Dim = 30;
  constexpr int Pop_Size = 100;
  constexpr int Pop_Size1 = 40;
  constexpr int Max = Dim * 10000;
  constexpr int Runs = 30;
  char dir[] = "./data";
  using namespace sevobench;
  evo_bench<CMA_ES, CEC2020, false, Dim, Pop_Size, Max, Runs, true, double>()
      .write_table(dir);
  evo_bench<DE, CEC2020, false, Dim, Pop_Size, Max, Runs, true, double>()
      .write_table(dir);
  evo_bench<SPSO2011, CEC2020, false, Dim, Pop_Size1, Max, Runs, true, double>()
      .write_table(dir);
  return 0;
}
