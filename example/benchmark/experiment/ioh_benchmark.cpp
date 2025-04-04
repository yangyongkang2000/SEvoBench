#include "de.h"
#include "ioh.hpp"
#include "mean_std.h"
#include <cstdio>
int main() {
  constexpr int Dim = 20;
  constexpr int Size = 12;
  constexpr int Runs = 30;
  constexpr int MaxFES = 10000 * Dim;
  constexpr int Step = 10 * Dim;
  std::vector<int> prob_idx(Size);
  std::iota(prob_idx.begin(), prob_idx.end(), 1001);
  auto suite = std::make_shared<ioh::suite::Suite<ioh::problem::CEC2022>>(
      prob_idx, std::vector<int>{1}, std::vector<int>{Dim}, "CEC2022");
  auto logger = std::make_shared<ioh::logger::Store>(
      std::vector<std::reference_wrapper<ioh::logger::Trigger>>{
          ioh::trigger::each(Step)},
      std::vector<std::reference_wrapper<ioh::logger::Property>>{
          ioh::watch::evaluations, ioh::watch::transformed_y_best});
  ioh::Experimenter<ioh::problem::CEC2022> experiment(
      suite, logger,
      [](const std::shared_ptr<ioh::problem::CEC2022> p) {
        de(*p, Dim, double(-100), double(100), MaxFES, 100, 0.5, 0.9);
      },
      Runs);
  auto t0 = std::chrono::high_resolution_clock ::now();
  experiment.run();
  auto t1 = std::chrono::high_resolution_clock ::now();
  for (int i = 0; i < Size; i++) {
    std::array<double, Runs> tmp;
    for (int j = 0; j < Runs; j++) {
      tmp[j] = logger
                   ->at(ioh::logger::Store::Cursor("CEC2022", 1001 + i, Dim, 1,
                                                   j, MaxFES / Step - 1),
                        "transformed_y_best")
                   .value();
    }
    auto [m, s] = mean_std(tmp.begin(), tmp.end());
    std::printf("F%d,mean:%f,std:%f\n", i + 1, m, s);
  }
  std::cout << "total time:"
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << "s\n";
  return 0;
}
