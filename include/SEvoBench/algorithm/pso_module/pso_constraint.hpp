#pragma once
#include "../../common/population.hpp"
#include "pso_velocity.hpp"
#include <algorithm>
#include <memory>

namespace sevobench::pso_module {

template <std::floating_point T> class velocity_constraint {
public:
  virtual void repair_velocity(particle_velocity<T> &) {}
  virtual ~velocity_constraint() {}
};

template <std::floating_point T> class positions_constraint {
public:
  virtual void repair_position(solution<T> &sol, T lb, T ub) {
    for (auto &_ : sol)
      _ = std::clamp(_, lb, ub);
  }
  virtual ~positions_constraint() = default;
};

template <std::floating_point T> struct pso_constraint {
  std::unique_ptr<velocity_constraint<T>> vc =
      std::make_unique<velocity_constraint<T>>();
  std::unique_ptr<positions_constraint<T>> pc =
      std::make_unique<positions_constraint<T>>();
  void repair_velocity(particle_velocity<T> &v) { vc->repair_velocity(v); }
  void repair_position(solution<T> &sol, T lb, T ub) {
    pc->repair_position(sol, lb, ub);
  }
};

template <std::floating_point T>
class spso_velocity_constraint : public velocity_constraint<T> {
  const T vmin;
  const T vmax;

public:
  spso_velocity_constraint(T _vmin, T _vmax) : vmin(_vmin), vmax(_vmax) {}
  void repair_velocity(particle_velocity<T> &v) override {
    for (auto &_ : v)
      _ = std::clamp(_, vmin, vmax);
  }
};

} // namespace sevobench::pso_module
