#include "SEvoBench/sevobench.hpp"
#include "de.h"
#include <iostream>
int main() {
  constexpr int Dim = 20;
  constexpr int MaxFES = 10000 * Dim;
  constexpr int Runs = 30;
  constexpr int Step = 10 * Dim;
  using T = double;
  auto suite = sevobench::problem::suite_builder<sevobench::problem::cec2022>()
                   .type<T>()
                   .dir(DATA_DIR "/cec2022_data/")
                   .dim<Dim>()
                   .problem_index(sevobench::problem::problem_range<1, 12>())
                   .build();
  sevobench::experiment::best_so_far_record<T> logger(suite, MaxFES, Runs,
                                                      Step);
  auto t0 = std::chrono::steady_clock::now();
  sevobench::experiment::evo_bench<PARALLEL>(
      [](auto &&f) {
        constexpr int PopSize = 100;
        de(f, Dim, f.lower_bound(), f.upper_bound(), MaxFES, PopSize,
           double(0.5), double(0.9));
      },
      suite, logger, Runs);
  auto t1 = std::chrono::steady_clock::now();
  auto m = logger.best().begin();
  for (auto i : suite.problem_index())
    for (int j = 0; j < suite.instance_count(); j++) {
      std::vector<T> tmp(m, m + Runs);
      m += Runs;
      auto s = sevobench::tool::mean_std(tmp.begin(), tmp.end());
      std::printf("F%d,instance:%d,mean:%f,std:%f,best:%f\n", i, j + 1, s[0],
                  s[1], s[2]);
    }
  std::cout << "total time:"
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << "s\n";
}