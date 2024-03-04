#include "SEvoBench/parallel_algorithm_benchmark.hpp"
#include <cstdio>
int main() {
  auto data = sevobench::evo_bench<sevobench::DE, sevobench::CEC2020>();
  auto table = data.get_table_data();
  auto alg_name = data.algorithm_name();
  auto pro_name = data.problem_name();
  std::printf("Algorithm:%s Suite:%s,Dim:%d\n", alg_name, pro_name,
              data.problem_dim());
  for (int i = 0; i < data.problem_size(); i++) {
    std::printf("F%d "
                "mean:%f,std:%f,best:%f,time:%fms,1/4best:%f,middle:%f,3/"
                "4best:%f,worst:%f\n",
                i + 1, table[8 * i], table[8 * i + 1], table[8 * i + 2],
                table[8 * i + 3], table[8 * i + 4], table[8 * i + 5],
                table[8 * i + 6], table[8 * i + 7]);
  }
  return 0;
}
