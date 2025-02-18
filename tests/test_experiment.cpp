#include "SEvoBench/sevobench.hpp"
#include <iostream>
int main() {
  constexpr int Dim = 30;
  constexpr int MaxFES = 10000 * Dim;
  constexpr int Runs = 10;
  using T = float;
  using namespace sevobench;
  auto suite = problem::suite_builder<problem::cec2017>()
                   .type<T>()
                   .instance_count(5)
                   .dim<Dim>()
                   .problem_index(problem::problem_range<20, 30>())
                   .build();
  sevobench::experiment::best_so_far_record<T> logger(suite, MaxFES, Runs,
                                                      100 * Dim);
  auto t0 = std::chrono::steady_clock::now();
  sevobench::experiment::evo_bench<true>(
      [](auto &&f) {
        using namespace sevobench::de_module;
        evolutionary_algorithm alg(MaxFES);
        population<T> pop(18 * Dim, Dim, f.lower_bound(), f.upper_bound());
        auto de =
            de_algorithm_builder<T>()
                .mutation(std::make_unique<ttpb1_weight_mutation<T>>(alg, 0.125,
                                                                     0.25))
                .parameter(std::make_unique<jso_parameter<T>>(alg))
                .crossover(std::make_unique<binomial_crossover<T>>())
                .constraint_handler(
                    std::make_unique<midpoint_target_repair<T>>())
                .population_strategy(
                    std::make_unique<linear_reduction<T>>(alg, 4, 18 * Dim))
                .archive(std::make_unique<fifo_archive<T>>())
                .build();
        de.run(pop, f, alg);
      },
      suite, logger, Runs);
  auto t1 = std::chrono::steady_clock::now();
  auto m = logger.best().begin();
  for (auto i : suite.problem_index())
    for (int j = 0; j < suite.instance_count(); j++) {
      {
        auto data = logger.at(i);
        auto data1 = logger.at(i, j + 1);
        for (int _ = 0; _ < Runs; _++) {
          if (std::abs(data[j][_].back().second - data1[_].back().second) >
              T(1e-3)) {
            std::printf("FAILED!\n");
            return -1;
          }
        }
      }
      std::vector<T> tmp(m, m + Runs);
      m += Runs;
      auto s = tool::mean_std(tmp.begin(), tmp.end());
      std::printf("F%d,instance:%d,mean:%f,std:%f,best:%f\n", i, j + 1,
                  s[0] - 100 * i, s[1], s[2] - 100 * i);
    }
  std::cout << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << "s\n";
}