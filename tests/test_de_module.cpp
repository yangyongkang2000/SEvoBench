#include "SEvoBench/sevobench.hpp"
#include <cstdio>
#include <iostream>
#include <map>
#include <string>

template <int N, typename T>
inline auto real_func(std::span<const T> t) noexcept {
  T sum = 0;
  for (int i = 0; i < N - 1; i++)
    sum += sevobench::tool::Pow<2>(t[i] - 1) +
           100 * sevobench::tool::Pow<2>((t[i + 1] - t[i] * t[i]));
  return sum;
}

template <std::floating_point T> inline auto generate_mutation() noexcept {
  using namespace sevobench::de_module;
  std::map<std::string, std::function<std::unique_ptr<de_mutation<T>>()>> table;
  table["rand1"] = [] { return std::make_unique<rand1_mutation<T>>(); };
  table["rand2"] = [] { return std::make_unique<rand2_mutation<T>>(); };
  table["best1"] = [] { return std::make_unique<best1_mutation<T>>(); };
  table["best2"] = [] { return std::make_unique<best2_mutation<T>>(); };
  table["ttpb1"] = [] { return std::make_unique<ttpb1_mutation<T>>(); };
  table["ttb1"] = [] { return std::make_unique<ttb1_mutation<T>>(); };
  table["ttb2"] = [] { return std::make_unique<ttb2_mutation<T>>(); };
  table["two_opt1"] = [] { return std::make_unique<two_opt1_mutation<T>>(); };
  table["two_opt2"] = [] { return std::make_unique<two_opt2_mutation<T>>(); };
  return table;
}

template <std::floating_point T> inline auto generate_crossover() noexcept {
  using namespace sevobench::de_module;
  std::map<std::string, std::function<std::unique_ptr<de_crossover<T>>()>>
      table;
  table["bin"] = [] { return std::make_unique<binomial_crossover<T>>(); };
  table["exp"] = [] { return std::make_unique<exponential_crossover<T>>(); };
  return table;
}

template <std::floating_point T> inline auto generate_constriant() noexcept {
  using namespace sevobench::de_module;
  std::map<std::string, std::function<std::unique_ptr<de_constraint<T>>()>>
      table;
  table["reflection"] = [] { return std::make_unique<reflection_repair<T>>(); };
  table["reinitial"] = [] {
    return std::make_unique<reinitialization_repair<T>>();
  };
  table["projection"] = [] { return std::make_unique<projection_repair<T>>(); };
  table["midpoint_target"] = [] {
    return std::make_unique<midpoint_target_repair<T>>();
  };
  table["midpoint_base"] = [] {
    return std::make_unique<midpoint_base_repair<T>>();
  };
  table["rand_base"] = [] { return std::make_unique<rand_base_repair<T>>(); };
  table["resample"] = [] { return std::make_unique<resample_repair<T>>(); };
  return table;
}

template <std::floating_point T>
inline auto generate_population_strategy() noexcept {
  using namespace sevobench::de_module;
  std::map<std::string, std::function<std::unique_ptr<de_population<T>>()>>
      table;
  table["default"] = [] { return std::make_unique<de_population<T>>(); };
  return table;
}

template <std::floating_point T> inline auto generate_parameter() noexcept {
  using namespace sevobench::de_module;
  std::map<std::string, std::function<std::unique_ptr<de_parameter<T>>()>>
      table;
  table["constant"] = [] { return std::make_unique<constant_parameter<T>>(); };
  table["jde"] = [] { return std::make_unique<jde_parameter<T>>(); };
  table["jade"] = [] { return std::make_unique<jade_parameter<T>>(); };
  table["shade"] = [] { return std::make_unique<shade_parameter<T>>(); };
  return table;
}

template <std::floating_point T> inline auto generate_art_parameter() noexcept {
  using namespace sevobench::de_module;
  std::map<std::string, std::function<std::unique_ptr<de_parameter<T>>()>>
      table;
  table["constant"] = [] { return std::make_unique<constant_parameter<T>>(); };
  table["jde"] = [] { return std::make_unique<jde_parameter<T>>(); };
  table["jade"] = [] { return std::make_unique<jade_parameter<T>>(); };
  table["shade"] = [] { return std::make_unique<shade_parameter<T>>(); };
  return table;
}

template <std::floating_point T> void test_basic_de() {
  using namespace sevobench::de_module;
  auto parameter_table = generate_parameter<T>();
  auto mutation_table = generate_mutation<T>();
  auto constraint_table = generate_constriant<T>();
  auto crossover_table = generate_crossover<T>();
  auto strategy_table = generate_population_strategy<T>();
  for (auto &[p_name, p] : parameter_table)
    for (auto &[m_name, m] : mutation_table)
      for (auto &[h_name, h] : constraint_table)
        for (auto &[c_name, c] : crossover_table)
          for (auto &[s_name, s] : strategy_table) {
            constexpr int Pop_Size = 100;
            constexpr int Dim = 30;
            sevobench::evolutionary_algorithm alg(Dim * 10000);
            sevobench::population<T> pop(Pop_Size, Dim, T(-100), T(100));
            auto de = de_algorithm_builder()
                          .mutation(m())
                          .parameter(p())
                          .constraint_handler(h())
                          .population_strategy(s())
                          .crossover(c())
                          .build();
            de.run(pop, real_func<Dim, T>, T(-100), T(100), alg);
            std::printf("%s_%s_%s_%s_%s"
                        ":",
                        p_name.c_str(), m_name.c_str(), h_name.c_str(),
                        c_name.c_str(), s_name.c_str());
            std::cout << std::min_element(pop.begin(), pop.end(),
                                          [](auto &x, auto &y) {
                                            return x.fitness() < y.fitness();
                                          })
                             ->fitness()
                      << '\n';
          }
}

template <std::floating_point T> void test_art_de() {
  using namespace sevobench::de_module;
  constexpr int Pop_Size = 671;
  constexpr int Dim = 30;
  sevobench::evolutionary_algorithm alg(Dim * 10000);
  sevobench::population<T> pop(Pop_Size, Dim, T(-100), T(100));
  auto de =
      de_algorithm_builder()
          .mutation(
              std::make_unique<ttpb1_weight_mutation<T>>(alg, 0.125, 0.25))
          .parameter(std::make_unique<jso_parameter<T>>(alg))
          .constraint_handler(std::make_unique<midpoint_target_repair<T>>())
          .crossover(std::make_unique<binomial_crossover<T>>())
          .archive(std::make_unique<fifo_archive<T>>());
  de.population_strategy(
        std::make_unique<non_linear_reduction<T>>(alg, 4, Pop_Size))
      .build()
      .run(pop, real_func<Dim, T>, T(-100), T(100), alg);
  std::cout << std::min_element(
                   pop.begin(), pop.end(),
                   [](auto &x, auto &y) { return x.fitness() < y.fitness(); })
                   ->fitness()
            << '\n';
}

int main() {
  test_basic_de<float>();
  test_art_de<float>();
  return 0;
}
