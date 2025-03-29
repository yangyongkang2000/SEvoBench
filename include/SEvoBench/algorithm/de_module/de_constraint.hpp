
#pragma once

#include "../../common/population.hpp"
#include "../../common/tool.hpp"
#include <memory>

namespace sevobench::de_module {
template <std::floating_point T> class de_constraint {
public:
  virtual void repair(std::span<T>, T, T, std::span<const T>,
                      std::span<const T>) {}

  virtual bool is_repair(std::span<T>, T, T, std::span<const T>,
                         std::span<const T>) {
    return true;
  }

  virtual ~de_constraint() = default;
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class reinitialization_repair final : public de_constraint<T> {
public:
  R RNG;

  void repair(std::span<T> p, T lb, T ub, std::span<const T>,
              std::span<const T>) override {
    for (auto &x : p)
      if (x < lb || x > ub)
        x = RNG.template rand_float<T>(lb, ub);
  }
};

template <std::floating_point T>
class projection_repair final : public de_constraint<T> {
public:
  void repair(std::span<T> p, T lb, T ub, std::span<const T>,
              std::span<const T>) override {
    for (auto &x : p)
      x = std::clamp(x, lb, ub);
  }
};

template <std::floating_point T>
class reflection_repair final : public de_constraint<T> {
public:
  void repair(std::span<T> p, T lb, T ub, std::span<const T>,
              std::span<const T>) override {
    for (auto &x : p)
      while (x < lb || x > ub) {
        if (x < lb) {
          x = 2 * lb - x;
        } else {
          x = 2 * ub - x;
        }
      }
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class rand_base_repair final : public de_constraint<T> {
public:
  R RNG;

  void repair(std::span<T> p, T lb, T ub, std::span<const T> base,
              std::span<const T>) override {
    for (std::size_t i = 0; i < p.size(); i++) {
      if (p[i] < lb) {
        p[i] = RNG.rand_float(lb, base[i]);
      } else if (p[i] > ub) {
        p[i] = RNG.rand_float(base[i], ub);
      }
    }
  }
};

template <std::floating_point T>
class midpoint_base_repair final : public de_constraint<T> {
public:
  void repair(std::span<T> p, T lb, T ub, std::span<const T> base,
              std::span<const T>) override {
    for (std::size_t i = 0; i < p.size(); i++) {
      if (p[i] < lb) {
        p[i] = T(0.5) * (lb + base[i]);
      } else if (p[i] > ub) {
        p[i] = T(0.5) * (ub + base[i]);
        ;
      }
    }
  }
};

template <std::floating_point T>
class midpoint_target_repair final : public de_constraint<T> {
public:
  void repair(std::span<T> p, T lb, T ub, std::span<const T>,
              std::span<const T> target) override {
    for (std::size_t i = 0; i < p.size(); i++) {
      if (p[i] < lb) {
        p[i] = T(0.5) * (lb + target[i]);
      } else if (p[i] > ub) {
        p[i] = T(0.5) * (ub + target[i]);
        ;
      }
    }
  }
};

template <std::floating_point T>
class resample_repair final : public de_constraint<T> {
  const int max_resamples = 100;
  std::unique_ptr<de_constraint<T>> ch =
      std::make_unique<projection_repair<T>>();
  int resamples = 0;

public:
  resample_repair() = default;

  resample_repair(std::unique_ptr<de_constraint<T>> &&_ch,
                  int _max_resamples = 100)
      : ch(std::move(_ch)), max_resamples(_max_resamples) {}

  bool is_repair(std::span<T> p, T lb, T ub, std::span<const T> base,
                 std::span<const T> target) override {
    if (resamples >= max_resamples) {
      resamples = 0;
      ch->repair(p, lb, ub, base, target);
      return true;
    }
    for (auto &x : p)
      if (x < lb || x > ub) {
        resamples++;
        return false;
      }
    return true;
  }
};

} // namespace sevobench::de_module
