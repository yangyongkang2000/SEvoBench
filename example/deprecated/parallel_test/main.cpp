#include "SEvoBench/parallel_algorithm_benchmark.hpp"
#include <array>
#include <chrono>
int main() {
  using namespace sevobench;
  constexpr int Dim = 10;
  constexpr int Pop_Size = 18 * Dim;
  constexpr int Max = Dim * 1000;
  constexpr int Runs = 30;
  constexpr int Repeat = 10;
  std::array<double, Repeat> t;
  for (int i = 0; i < Repeat; i++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    evo_bench<LSHADE, CEC2022, false, Dim, Pop_Size, Max, Runs, PARALLEL>();
    auto t1 = std::chrono::high_resolution_clock::now();
    t[i] =
        (double)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
            .count();
  }
  std::ofstream os(PARALLEL ? "./parallel.txt" : "./single.txt");
  for (auto _ : t)
    os << _ << '\n';
}
