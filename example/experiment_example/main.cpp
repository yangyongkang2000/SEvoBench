#include <SEvoBench/sevobench.hpp>

int main() {
  using namespace sevobench;
  using namespace experiment;
  using namespace de_module;
  // 1. Configure test suite
  auto suite = problem::suite_builder<sevobench::problem::cec2017>()
                   .template dim<30>()
                   .type<double>()
                   .problem_index(problem::problem_range<1, 15>())
                   .instance_count(5)
                   .build();
  // 2. Initialize performance tracker
  best_so_far_record<double> recorder(suite, 1e5, 25, 1000);

  // 4. Execute benchmark
  evo_bench(
      [](auto &p) {
        evolutionary_algorithm tracker(1e5, 100, 30);
        // 3. Configure JADE algorithm
        auto config = de_config<true, double>{
            std::make_unique<jade_parameter<double>>(0.1, 0.5),
            std::make_unique<ttpb1_mutation<double>>(0.11),
            std::make_unique<projection_repair<double>>(),
            std::make_unique<binomial_crossover<double>>(),
            std::make_unique<linear_reduction<double>>(tracker, 50, 100),
            std::make_unique<fifo_archive<double>>(2.0)};
        population<double> pop(100, tracker.dim(), double(-100), double(100));
        auto de = de_algorithm(std::move(config));
        return de.run(
            pop, [&](std::span<const double> x) { return p(x); }, double(-100),
            double(100), tracker);
      },
      suite, recorder,
      25 // 25 independent runs
  );

  // 5. Export results
  std::ofstream csv("results.csv");
  csv << "Problem,Instance,Run,FEs,Fitness\n";
  for (size_t pid : suite.problem_index()) {
    for (int inst = 1; inst <= suite.instance_count(); inst++) {
      auto data = recorder.at(pid, inst);
      for (size_t run = 0; run < data.size(); run++) {
        for (const auto &[fes, fit] : data[run]) {
          csv << pid << "," << inst << "," << run + 1 << "," << fes << ","
              << fit << "\n";
        }
      }
    }
  }
}