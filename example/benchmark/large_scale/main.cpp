#include "SEvoBench/problem/large_scale/cec2010.hpp"
#include "SEvoBench/sevobench.hpp"
#include <cassert>
#include <chrono>
#include <iostream>
int main() {
  constexpr int Pop_Size = 300;
  constexpr int Group = 50;
  constexpr int Dim = 1000;
  constexpr int MaxFES = 1000 * Dim;
  constexpr int Runs = 5;
  constexpr int Step = 1000 * Dim;
  using T = float;
  using namespace sevobench;
  auto suite = problem::suite_builder<
                   problem::simd::cec2010_setting<Group, USING_SIMD>::cec2010>()
                   .type<T>()
                   .instance_count(1)
                   .dim<Dim>()
                   .problem_index(Full ? problem::problem_range<1, 20>()
                                       : problem::problem_range<1, 6>())
                   .build();
  sevobench::experiment::best_so_far_record<T> logger(suite, MaxFES, Runs,
                                                      Step);
  auto t0 = std::chrono::steady_clock::now();
  sevobench::experiment::evo_bench<Parallel>(
      [](auto &&f) {
        other_algorithm::cso<Dim, Pop_Size, MaxFES>(
            [&](const T *x) { return f(std::span<const T>(x, Dim)); },
            f.lower_bound(), f.upper_bound());
      },
      suite, logger, Runs);
  auto t1 = std::chrono::steady_clock::now();
  auto m = logger.best().begin();
  for (auto i : suite.problem_index())
    for (int j = 0; j < suite.instance_count(); j++) {
      std::vector<T> tmp(m, m + Runs);
      m += Runs;
      auto s = tool::mean_std(tmp.begin(), tmp.end());
      std::printf("F%d,instance:%d,mean:%f,std:%f,best:%f\n", i, j + 1, s[0],
                  s[1], s[2]);
    }
  std::cout << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << "s\n";
  return 0;
}
