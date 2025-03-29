#pragma once
#include "tool.hpp"
#include <algorithm>
#include <concepts>
#include <vector>
namespace sevobench {

template <std::floating_point T = float>
class solution : public std::vector<T> {
  T _fitness = std::numeric_limits<T>::max();

public:
  using std::vector<T>::vector;
  auto fitness() const noexcept { return _fitness; }
  auto set_fitness(T f) noexcept { _fitness = f; }
  auto dim() const noexcept { return static_cast<int>(this->size()); }
  auto evaluate(auto &&f) noexcept {
    constexpr auto b1 = requires {
      { f(*this) } -> std::same_as<T>;
    };
    constexpr auto b2 = requires {
      { f(this->data()) } -> std::same_as<T>;
    };
    constexpr auto b3 = requires {
      { f(this->data(), this->size()) } -> std::same_as<T>;
    };
    if constexpr (requires {
                    { f(*this) } -> std::same_as<T>;
                  }) {
      _fitness = f(*this);
    } else if constexpr (requires {
                           { f(this->data()) } -> std::same_as<T>;
                         }) {
      _fitness = f(this->data());
    } else if constexpr (requires {
                           { f(this->data(), this->size()) } -> std::same_as<T>;
                         }) {
      _fitness = f(this->data(), this->size());
    } else {
      static_assert(b1 || b2 || b3, "TYPE OF f IS INVALID");
    }
    return _fitness;
  }
  template <typename R = tool::rng>
    requires tool::random_generator_concept<R, T>
  auto randomize(T lb, T ub, R &&RNG = R()) noexcept {
    std::generate(this->begin(), this->end(),
                  [=, &RNG] { return RNG.rand_float(lb, ub); });
  }
};

template <std::floating_point T = float>
class population : public std::vector<solution<T>> {
  int _dim = 0;

public:
  template <typename R = tool::rng>
    requires tool::random_generator_concept<R, T>
  auto randomize(T lb, T ub, R &&RNG = R()) noexcept {
    for (auto &_ : *this) {
      _.randomize(lb, ub, RNG);
    }
  }
  population() = default;
  population(int dim_) : _dim(dim_) {}
  population(int _pop_size, int dim_)
      : std::vector<solution<T>>(_pop_size, solution<T>(dim_)), _dim(dim_) {}

  template <typename R = tool::rng>
    requires tool::random_generator_concept<R, T>
  population(int s, int d, T lb, T ub, R &&RNG = R()) : population<T>(s, d) {
    randomize(lb, ub, RNG);
  }
  auto pop_size() const noexcept { return static_cast<int>(this->size()); }
  auto set_dim(int dim_) noexcept { _dim = dim_; }
  auto dim() const noexcept { return _dim; }
};
} // namespace sevobench
