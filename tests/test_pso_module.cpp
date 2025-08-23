#include "SEvoBench/sevobench.hpp"
#include <cstdio>
#include <iostream>
#include <map>
#include <string>

template <int N, typename T>
inline auto real_func(std::span<const T> t)

    noexcept {
  T sum = 0;
  for (int i = 0; i < N - 1; i++)
    sum += sevobench::tool::Pow<2>(t[i] - 1) +
           100 * sevobench::tool::Pow<2>((t[i + 1] - t[i] * t[i]));
  return sum;
}

template <std::floating_point T>
inline auto generate_update()

    noexcept {
  using namespace sevobench::pso_module;
  std::map<std::string, std::function<std::unique_ptr<pso_update<T>>()>> table;
  table["inertia_weight"] = [] {
    return

        std::make_unique<inertia_weight_update<T>>();
  };
  table["bare_bone"] = [] {
    return

        std::make_unique<bare_bone_update<T>>();
  };
  table["fips"] = [] {
    return

        std::make_unique<fips_update<T>>();
  };
  table["spherical"] = [] {
    return

        std::make_unique<spherical_update<T>>();
  };
  return table;
}

template <std::floating_point T>
inline auto generate_topology()

    noexcept {
  using namespace sevobench::pso_module;
  std::map<std::string, std::function<std::unique_ptr<pso_topology<T>>()>>
      table;
  table["gbest"] = [] {
    return

        std::make_unique<gbest_topology<T>>();
  };
  table["lbest"] = [] {
    return

        std::make_unique<lbest_topology<T>>();
  };
  table["random"] = [] {
    return

        std::make_unique<random_topology<T>>();
  };
  table["von_neumann"] = [] {
    return

        std::make_unique<von_neumann_topology<T>>();
  };
  table["increasing"] = [] {
    return std::make_unique<increasing_topology<T>>(277);
  };
  table["decreasing"] = [] {
    return std::make_unique<decreasing_topology<T>>(277);
  };
  table["dms"] = [] { return std::make_unique<dms_topology<T>>(5); };
  return table;
}

template <std::floating_point T> void test_pso() {
  for (auto &[u_name, u] : generate_update<T>())
    for (auto &[t_name, t] : generate_topology<T>()) {
      constexpr int Pop_Size = 40;
      constexpr int Dim = 30;
      sevobench::evolutionary_algorithm alg(10000 * Dim);
      sevobench::population<T> pop(Pop_Size, Dim, T(-100), T(100));
      sevobench::pso_module::pso_velocity<T> vec(
          Pop_Size, sevobench::pso_module::particle_velocity<T>(Dim));
      auto pso =
          sevobench::pso_module::pso_algorithm_builder<T>()
              .update(u())
              .topology(t())
              .constraint_handler(sevobench::pso_module::pso_constraint<T>{
                  .vc = std::make_unique<
                      sevobench::pso_module::spso_velocity_constraint<T>>(-40,
                                                                          40)})
              .build();
      pso.replace_update(u());
      pso.replace_topology(t());
      pso.replace_constraint_handler(sevobench::pso_module::pso_constraint<T>{
          .vc = std::make_unique<
              sevobench::pso_module::spso_velocity_constraint<T>>(-40, 40)});
      pso.run(pop, vec, real_func<Dim, T>, T(-100), T(100), alg);
      std::printf("%s_%s:%f\n", u_name.c_str(), t_name.c_str(),
                  pso.topology()->best_value());
    }
}

int main() { test_pso<float>(); }