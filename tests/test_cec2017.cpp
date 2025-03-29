#include <iostream>

#include "SEvoBench/sevobench.hpp"

void cec17_test_func(double *, double *, int, int, int);

double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag, *SS;

int main() {
  using namespace sevobench;
  constexpr int Dim = 30;
  std::array<double, Dim> time{};
  std::array<double, Dim> time1{};
  std::array<double, Dim> f1{};
  std::array<double, Dim> f2{};
  population<double> pop(30, Dim, double(-100), double(100));
  {
    tool::simple_rand sr{100};
    auto suites = problem::suite_builder<problem::cec2017>()
                      .type<double>()
                      .dir(DATA_DIR "/cec2017_data/")
                      .dim<Dim>()
                      .problem_index(problem::problem_range<1, 30>())
                      .build();
    int i = 0;
    for (auto &p : suites) {
      auto t1 = std::chrono::high_resolution_clock::now();
      for (int _ = 0; _ < 1000; _++)
        f1[i] = (*p)(pop[sr() % 30]);
      auto t2 = std::chrono::high_resolution_clock::now();
      time[i++] = (t2 - t1).count();
    }
  }
  {
    tool::simple_rand sr{100};
    for (int i = 0; i < 30; i++) {
      auto t1 = std::chrono::high_resolution_clock::now();
      for (int _ = 0; _ < 1000; _++)
        cec17_test_func(pop[sr() % 30].data(), &f2[i], Dim, 1, i + 1);
      auto t2 = std::chrono::high_resolution_clock::now();
      time1[i] = (t2 - t1).count();
    }
  }
  for (int i = 0; i < 30; i++) {
    printf("%d,delta:%f,ratio:%f\n", i + 1, std::abs(f1[i] - f2[i]),
           time1[i] / time[i]);
  }
  std::cout << "ratio:"
            << std::accumulate(time1.begin(), time1.end(), double(0)) /
                   std::accumulate(time.begin(), time.end(), double(0))
            << "\n";
  return 0;
}
