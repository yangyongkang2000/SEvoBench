#include "SEvoBench/parallel_algorithm_benchmark.hpp"
#include <cassert>
template <std::uint64_t Alg, std::uint64_t Prob> void test_benchmark() {
  auto r1 = sevobench::evo_bench<Alg, Prob>();
  auto m = r1.get_table_data();
  assert(m.size() == 8 * r1.size);
  for (int i = 0; i < r1.size; i++) {
    assert(m[8 * i + 2] <= m[8 * i + 4]);
    for (int j = 4; j < 6; j++)
      assert(m[8 * i + j] <= m[8 * i + j + 1]);
  }
  auto r2 = sevobench::evo_bench<Alg, Prob, true>();
  m = r2.get_table_data();
  assert(m.size() == 8 * r2.size);
  for (int i = 0; i < r2.size; i++) {
    assert(m[8 * i + 2] <= m[8 * i + 4]);
    for (int j = 4; j < 6; j++)
      assert(m[8 * i + j] <= m[8 * i + j + 1]);
  }
  auto w1 = r2.get_curve_data();
  assert(w1.size() == 2 * r2.max * r2.size);
  for (int i = 0; i < r2.size; i++) {
    for (int j = 0; j < r2.max - 1; j++) {
      assert(w1[2 * r2.max * i + 2 * j] < w1[2 * r2.max * i + 2 * (j + 1)]);
      assert(w1[2 * r2.max * i + 2 * j + 1] >=
             w1[2 * r2.max * i + 2 * (j + 1) + 1]);
    }
  }
  r2.write_curve("./curve_data");
  r1.write_table("./table_data");
}
int main() {
  test_benchmark<sevobench::PSO, sevobench::ShiftFunc>();
  test_benchmark<sevobench::SLPSO, sevobench::ShiftFunc>();
  test_benchmark<sevobench::CSO, sevobench::ShiftFunc>();
  test_benchmark<sevobench::SPSO2007, sevobench::ShiftFunc>();
  test_benchmark<sevobench::SPSO2011, sevobench::ShiftFunc>();
  test_benchmark<sevobench::DE, sevobench::ShiftFunc>();
  test_benchmark<sevobench::JADE, sevobench::ShiftFunc>();
  test_benchmark<sevobench::SHADE, sevobench::ShiftFunc>();
  test_benchmark<sevobench::LSHADE, sevobench::ShiftFunc>();
  test_benchmark<sevobench::ABC, sevobench::ShiftFunc>();
  test_benchmark<sevobench::ES, sevobench::ShiftFunc>();
  test_benchmark<sevobench::RandomSearch, sevobench::ShiftFunc>();
}