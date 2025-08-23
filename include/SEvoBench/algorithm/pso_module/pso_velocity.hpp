#pragma once
#include <concepts>
#include <vector>
namespace sevobench::pso_module {
template <std::floating_point T> using particle_velocity = std::vector<T>;
template <std::floating_point T>
using pso_velocity = std::vector<particle_velocity<T>>;
} // namespace sevobench::pso_module
