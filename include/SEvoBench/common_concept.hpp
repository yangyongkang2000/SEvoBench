#pragma once
#include <concepts>
#include <type_traits>
namespace sevobench {
template <auto Dim, auto Pop_Size, auto Max>
concept algorithm_parameter_concept = Dim >= 1 && Pop_Size >= 1 && Max >= 1;

template <typename F, typename T>
concept algorithm_problem_concept = requires(F &&f, T *x) {
  requires std::same_as<T, std::remove_cvref_t<decltype(f(x))>>;
};

template <typename F, typename G, typename T>
concept algorithm_positions_concept = requires(G &&positions) {
  requires std::same_as<std::remove_reference_t<decltype(positions[0][0])>, T>;
  { positions.begin() };
  { positions.end() };
  requires algorithm_problem_concept<
      F, std::remove_pointer_t<decltype((*positions.begin()).data())>>;
  requires algorithm_problem_concept<
      F, std::remove_pointer_t<decltype(positions[0].data())>>;
};
template <typename F, typename G, typename T>
concept algorithm_vector_concept = requires(G &&positions) {
  requires std::same_as<std::remove_reference_t<decltype(positions[0])>, T>;
  requires algorithm_problem_concept<
      F, std::remove_pointer_t<decltype(positions.data())>>;
};

template <auto Dim, auto Pop_Size, auto Max, typename F, typename T>
concept algorithm_func_concept =
    algorithm_parameter_concept<Dim, Pop_Size, Max> &&
    algorithm_problem_concept<F, T>;

} // namespace sevobench
