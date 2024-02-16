#pragma once
#include "abc.hpp"
#include "cso.hpp"
#include "de.hpp"
#include "es.hpp"
#include "jade.hpp"
#include "pso.hpp"
#include "random_search.hpp"
#include "slpso.hpp"
#include "template_pso.hpp"
#include "template_shade.hpp"
#include "tool.hpp"

#define SEvoBench_REGISTER_ALGORITHM(func_name, func, parameter)               \
  namespace sevobench {                                                        \
  static constexpr auto func_name =                                            \
      tool::alg_hash<std::uint64_t>("sevobench" #func);                        \
  template <int Dim, int Pop_Size, int Max, bool Memory_Flag>                  \
    requires algorithm_parameter_concept<Dim, Pop_Size, Max>                   \
  struct single_algorithm<func_name, Dim, Pop_Size, Max, Memory_Flag> {        \
    static constexpr char name[] = #func_name;                                 \
    static constexpr int dim = Dim;                                            \
    static constexpr int size = Pop_Size;                                      \
    static constexpr int max = Max;                                            \
    static constexpr bool flag = Memory_Flag;                                  \
    template <typename F, std::floating_point T,                               \
              typename Parameter_Type = parameter<T>>                          \
      requires algorithm_problem_concept<F, T>                                 \
    inline auto                                                                \
    operator()(F &&f, T l, T r,                                                \
               const Parameter_Type &pt = Parameter_Type()) noexcept {         \
      return func<Dim, Pop_Size, Max, Memory_Flag>(f, l, r, pt);               \
    }                                                                          \
  };                                                                           \
  template <std::floating_point T>                                             \
  struct single_algorithm_parameter<func_name, T> : public parameter<T> {      \
    using base_type = parameter<T>;                                            \
  };                                                                           \
  }

namespace sevobench {
template <std::uint64_t Alg_HashName, int Dim, int Pop_Size, int Max,
          bool Memory_Flag>
  requires algorithm_parameter_concept<Dim, Pop_Size, Max>
struct single_algorithm;

template <std::uint64_t Alg_HashName, std::floating_point T>
struct single_algorithm_parameter;

} // namespace sevobench

SEvoBench_REGISTER_ALGORITHM(SLPSO, slpso, slpso_parameter)
    SEvoBench_REGISTER_ALGORITHM(RandomSearch, random_search,
                                 random_search_parameter)
        SEvoBench_REGISTER_ALGORITHM(PSO, pso, pso_parameter)
            SEvoBench_REGISTER_ALGORITHM(ABC, abc, abc_parameter)
                SEvoBench_REGISTER_ALGORITHM(CSO, cso, cso_parameter)
                    SEvoBench_REGISTER_ALGORITHM(ES, es, es_parameter)
                        SEvoBench_REGISTER_ALGORITHM(JADE, jade, jade_parameter)
                            SEvoBench_REGISTER_ALGORITHM(SHADE, shade,
                                                         shade_parameter)
                                SEvoBench_REGISTER_ALGORITHM(LSHADE, lshade,
                                                             lshade_parameter)
                                    SEvoBench_REGISTER_ALGORITHM(DE, de,
                                                                 de_parameter)
                                        SEvoBench_REGISTER_ALGORITHM(
                                            SPSO2007, spso2007,
                                            spso2007_parameter)
                                            SEvoBench_REGISTER_ALGORITHM(
                                                SPSO2011, spso2011,
                                                spso2011_parameter)

#define REGISTER_ALGORITHM(func_name, func, parameter)                         \
  static constexpr auto func_name =                                            \
      sevobench::tool::alg_hash<std::uint64_t>(#func);                         \
  template <int Dim, int Pop_Size, int Max, bool Memory_Flag>                  \
    requires sevobench::algorithm_parameter_concept<Dim, Pop_Size, Max>        \
  struct sevobench::single_algorithm<func_name, Dim, Pop_Size, Max,            \
                                     Memory_Flag> {                            \
    static constexpr char name[] = #func_name;                                 \
    static constexpr int dim = Dim;                                            \
    static constexpr int size = Pop_Size;                                      \
    static constexpr int max = Max;                                            \
    static constexpr bool flag = Memory_Flag;                                  \
    template <typename F, std::floating_point T,                               \
              typename Parameter_Type = parameter<T>>                          \
      requires sevobench::algorithm_problem_concept<F, T>                      \
    inline auto                                                                \
    operator()(F &&f, T l, T r,                                                \
               const Parameter_Type &pt = Parameter_Type()) noexcept {         \
      return func<Dim, Pop_Size, Max, Memory_Flag, F, T>(f, l, r, pt);         \
    }                                                                          \
  };                                                                           \
  template <std::floating_point T>                                             \
  struct sevobench::single_algorithm_parameter<func_name, T>                   \
      : public parameter<T> {                                                  \
    using base_type = parameter<T>;                                            \
  };
