

#pragma once

#include "de_config.hpp"

namespace sevobench::de_module {

template <bool Use_Archive, std::floating_point T> class de_algorithm {
  de_config<Use_Archive, T> config_;

public:
  using value_type = T;
  static constexpr auto use_archive() { return Use_Archive; }
  de_algorithm(de_algorithm &&) = default;

  de_algorithm &operator=(de_algorithm &&) = default;

  de_algorithm(de_config<Use_Archive, T> &&_config)
      : config_(std::move(_config)) {}

  template <typename F = evolutionary_algorithm_condition>
  void run(std::convertible_to<population<T>> auto &&pop, auto &&f, T lb, T ub,
           evolutionary_algorithm &alg, F condition = {})

      noexcept {
    population<T> trial(pop.pop_size(), pop.dim());
    alg.set_pop_size(pop.pop_size());
    alg.set_dim(pop.dim());
    if constexpr (Use_Archive) {
      config_.archive->prepare(pop);
    }
    for (auto &_ : pop)
      _.evaluate(f);
    alg.add_fes(pop.pop_size());
    do {
      config_.iterator(pop, trial, f, lb, ub, alg);
    } while (condition(pop, alg));
  }
  template <typename F = evolutionary_algorithm_condition>
  void run(std::convertible_to<population<T>> auto &&pop, auto &&f,
           evolutionary_algorithm &alg, F condition = {})

      noexcept {
    run(pop, f, f.lower_bound(), f.upper_bound(), alg, condition);
  }

  [[nodiscard]] auto mutation() noexcept { return config_.mutation.get(); }
  auto replace_mutation(std::unique_ptr<de_mutation<T>> &&ptr) noexcept {
    config_.mutation = std::move(ptr);
  }
  auto swap_mutation(std::unique_ptr<de_mutation<T>> &ptr) noexcept {
    config_.mutation.swap(ptr);
  }
  [[nodiscard]] auto parameter() noexcept { return config_.parameter.get(); }
  auto replace_parameter(std::unique_ptr<de_parameter<T>> &&ptr) noexcept {
    config_.parameter = std::move(ptr);
  }
  auto swap_parameter(std::unique_ptr<de_parameter<T>> &ptr) noexcept {
    config_.parameter.swap(ptr);
  }
  [[nodiscard]] auto constraint_handler() noexcept {
    return config_.constraint_handler.get();
  }
  auto
  replace_constraint_handler(std::unique_ptr<de_constraint<T>> &&ptr) noexcept {
    config_.constraint_handler = std::move(ptr);
  }
  auto
  swap_constraint_handler(std::unique_ptr<de_constraint<T>> &ptr) noexcept {
    config_.constraint_handler.swap(ptr);
  }
  [[nodiscard]] auto crossover() noexcept { return config_.crossover.get(); }
  auto replace_crossover(std::unique_ptr<de_crossover<T>> &&ptr) noexcept {
    config_.crossover = std::move(ptr);
  }
  auto swap_crossover(std::unique_ptr<de_crossover<T>> &ptr) noexcept {
    config_.crossover.swap(ptr);
  }
  [[nodiscard]] auto population_strategy() noexcept {
    return config_.population_strategy.get();
  }
  auto replace_population_strategy(
      std::unique_ptr<de_population<T>> &&ptr) noexcept {
    config_.population_strategy = std::move(ptr);
  }
  auto
  swap_population_strategy(std::unique_ptr<de_population<T>> &ptr) noexcept {
    config_.population_strategy.swap(ptr);
  }
  auto archive() noexcept {
    if constexpr (Use_Archive) {
      return config_.archive.get();
    }
  }
  auto replace_archive([[maybe_unused]] std::conditional_t<
                       Use_Archive, std::unique_ptr<de_archive<T>> &&, void *>
                           ptr) noexcept {
    if constexpr (Use_Archive) {
      config_.archive = std::move(ptr);
    }
  }
  auto swap_archive([[maybe_unused]] std::conditional_t<
                    Use_Archive, std::unique_ptr<de_archive<T>> &, void *>
                        ptr) noexcept {
    if constexpr (Use_Archive) {
      config_.archive.swap(ptr);
    }
  }
  void iterator(auto &&pop, auto &&trial, auto &&f, T lb, T ub,
                evolutionary_algorithm &alg) noexcept {
    config_.iterator(pop, trial, f, lb, ub, alg);
  }
};

template <std::floating_point T = float, bool Use_Archive = false,
          bool B1 = false, bool B2 = false, bool B3 = false, bool B4 = false,
          bool B5 = false>
class [[nodiscard]] de_algorithm_builder {
  de_config<Use_Archive, T> config_;

public:
  de_algorithm_builder() = default;

  de_algorithm_builder(de_config<Use_Archive, T> &&_config)
      : config_(std::move(_config)) {}

  [[nodiscard]] auto mutation(std::unique_ptr<de_mutation<T>> &&m)

      noexcept {
    config_.mutation = std::move(m);
    return de_algorithm_builder<T, Use_Archive, true, B2, B3, B4, B5>(
        std::move(config_));
  }

  [[nodiscard]] auto constraint_handler(std::unique_ptr<de_constraint<T>> &&c)

      noexcept {
    config_.constraint_handler = std::move(c);
    return de_algorithm_builder<T, Use_Archive, B1, true, B3, B4, B5>(
        std::move(config_));
  }

  [[nodiscard]] auto crossover(std::unique_ptr<de_crossover<T>> &&c)

      noexcept {
    config_.crossover = std::move(c);
    return de_algorithm_builder<T, Use_Archive, B1, B2, true, B4, B5>(
        std::move(config_));
  }

  [[nodiscard]] auto parameter(std::unique_ptr<de_parameter<T>> &&p)

      noexcept {
    config_.parameter = std::move(p);
    return de_algorithm_builder<T, Use_Archive, B1, B2, B3, true, B5>(
        std::move(config_));
  }

  [[nodiscard]] auto population_strategy(std::unique_ptr<de_population<T>> &&p)

      noexcept {
    config_.population_strategy = std::move(p);
    return de_algorithm_builder<T, Use_Archive, B1, B2, B3, B4, true>(
        std::move(config_));
  }

  [[nodiscard]] auto archive(std::unique_ptr<de_archive<T>> &&a)

      noexcept {
    de_config<true, T> _config{
        .parameter = std::move(config_.parameter),
        .mutation = std::move(config_.mutation),
        .constraint_handler = std::move(config_.constraint_handler),
        .crossover = std::move(config_.crossover),
        .population_strategy = std::move(config_.population_strategy),
        .archive = std::move(a)};
    return de_algorithm_builder<T, true, B1, B2, B3, B4, B5>(
        std::move(_config));
  }

  [[nodiscard]] auto build()

      noexcept {
    static_assert(B1, "MUTATION IS MISSING!");
    static_assert(B2, "CONSTRAINT HANDLER IS MISSING!");
    static_assert(B3, "CROSSOVER IS MISSING!");
    static_assert(B4, "MUTATION IS MISSING!");
    static_assert(B5, "POPULATION STRATEGY IS MISSING!");
    return de_algorithm(std::move(config_));
  }
};

} // namespace sevobench::de_module
