#pragma once
#include "parallel_task.hpp"
#include "single_algorithm.hpp"
#include "single_problem.hpp"
#include "tool.hpp"
#include <chrono>
#include <deque>
#include <filesystem>
#include <fstream>

namespace sevobench {

namespace sevobench_detail {

template <typename F, typename Algorithm, bool Memory_Flag = false,
          typename Parameter_Type>
inline auto
benchamrk_calc(const Parameter_Type &par = Parameter_Type()) noexcept {
  if constexpr (Memory_Flag) {
    auto [_, u, v] = Algorithm()(F(), F().L, F().U, par);
    return std::make_pair(u, std::move(v));
  }
  if constexpr (!Memory_Flag) {
    auto [_, v] = Algorithm()(F(), F().L, F().U, par);
    return v;
  }
}

template <typename V, typename T>
concept curve_concept = std::movable<V> && requires(V v) {
  requires std::same_as<T, std::remove_reference_t<decltype(v[0])>>;
};

template <std::uint64_t Alg_HashName, std::uint64_t Prob_HashName, int Dim,
          int Pop_Size, int Max, bool Memory_Flag, int Prob_Index,
          typename Parameter_Type, typename T,
          typename Alg =
              single_algorithm<Alg_HashName, Dim, Pop_Size, Max, Memory_Flag>,
          typename Prob = single_problem<Prob_HashName, Prob_Index, Dim, T>>
concept single_algorithm_problem_concept = requires(const Parameter_Type &pt) {
  requires algorithm_parameter_concept<Dim, Pop_Size, Max>;
  requires Prob_Index >= 1;
  requires Prob::dim == Dim &&Prob::index == Prob_Index &&Prob::size >= 1 &&
                   Prob::L < Prob::U;
  requires Alg::dim == Dim &&Alg::flag == Memory_Flag &&Alg::max == Max;
  Alg::name;
  Alg()(Prob(), Prob::L, Prob::U, pt);
  (std::tuple_size_v<decltype(Alg()(Prob(), Prob::L, Prob::U, pt))>) ==
      (Memory_Flag ? 3 : 2);
  requires std::same_as<
      T,
      std::tuple_element_t<1, decltype(Alg()(Prob(), Prob::L, Prob::U, pt))>>;
  requires curve_concept<
      std::tuple_element_t<2, decltype(single_algorithm<Alg_HashName, Dim,
                                                        Pop_Size, Max, true>()(
                                  Prob(), Prob::L, Prob::U, pt))>,
      T>;
  requires std::same_as<
      bool, std::remove_reference_t<
                decltype(init_single_problem<Prob_HashName, Dim, T>::is_init)>>;
  init_single_problem<Prob_HashName, Dim, T>{};
};

template <std::uint64_t Alg_HashName, std::uint64_t Prob_HashName, int Dim,
          int Pop_Size, int Max, bool Memory_Flag, int Prob_Index,
          typename Parameter_Type, typename T, typename Task,
          typename Alg =
              single_algorithm<Alg_HashName, Dim, Pop_Size, Max, Memory_Flag>,
          typename Prob = single_problem<Prob_HashName, Prob_Index, Dim, T>>
concept single_algorithm_problem_task_concept = requires(
    const Parameter_Type &par, Task &pool) {
  requires single_algorithm_problem_concept<
      Alg_HashName, Prob_HashName, Dim, Pop_Size, Max, Memory_Flag, Prob_Index,
      Parameter_Type, T, Alg, Prob>;
  {
    pool.submit(benchamrk_calc<Prob, Alg, Memory_Flag, Parameter_Type>, par)
  } -> std::same_as<std::future<std::conditional_t<
      Memory_Flag,
      std::pair<T, typename std::tuple_element_t<Memory_Flag ? 2 : 0,
                                                 decltype(Alg()(Prob(), Prob::L,
                                                                Prob::U))>>,
      T>>>;
};

template <std::floating_point T, typename Problem> class tmp_benchmark_result {
  static constexpr auto P = Problem::size;
  std::deque<T> graph_data;
  std::deque<T> table_data;
  template <bool b> auto _write(const char *file_name) {
    auto &m = b ? graph_data : table_data;
    if (!m.empty()) {
      std::ofstream os(file_name);
      os << P << '\n';
      for (auto x : m)
        os << x << '\n';
    }
  }

public:
  using value_type = T;
  auto write_table(const char *file_name) { _write<false>(file_name); }
  auto write_curve(const char *file_name) { _write<true>(file_name); }

  template <bool b, typename V> auto insert(const V &v) {
    if constexpr (b)
      table_data.insert(table_data.end(), std::begin(v), std::end(v));
    else
      graph_data.insert(graph_data.end(), std::begin(v), std::end(v));
  }
};

template <std::floating_point T, typename Algorithm, typename Problem>
class benchmark_result {
  std::vector<T> graph_data;
  std::vector<T> table_data;

public:
  using value_type = T;
  static constexpr auto algname = Algorithm::name;
  static constexpr auto proname = Problem::name;
  static constexpr auto pop = Algorithm::size;
  static constexpr auto dim = Algorithm::dim;
  static constexpr auto size = Problem::size;
  static constexpr auto flag = Algorithm::flag;
  static constexpr auto max = Algorithm::max;
  auto get_curve_data() {
    static_assert(flag, "Graph Data is Empty.");
    return graph_data;
  }
  auto get_table_data() { return table_data; }
  constexpr auto problem_size() { return size; }
  constexpr auto problem_name() { return proname; }
  constexpr auto algorithm_name() { return algname; }
  constexpr auto problem_dim() { return dim; }
  constexpr auto algorithm_max() { return max; }
  constexpr auto algorithm_flag() { return flag; }
  void write_table(const char *dir_name) {
    namespace fs = std::filesystem;
    fs::path dir(dir_name);
    if (!fs::is_directory(fs::status(dir))) {
      fs::create_directory(dir);
    }
    auto file_path =
        dir / (std::string(algname) + "_" + std::to_string(pop) + "_" +
               std::string(proname) + "_" + std::to_string(dim) + "_table.txt");
    if (!table_data.empty()) {
      std::ofstream os(file_path.string());
      for (auto x : table_data)
        os << x << '\n';
    }
  }
  void write_curve(const char *dir_name) {
    static_assert(flag, "Graph Data is Empty.");
    namespace fs = std::filesystem;
    fs::path dir(dir_name);
    if (!fs::is_directory(fs::status(dir))) {
      fs::create_directory(dir);
    }
    if (!graph_data.empty()) {
      auto len = graph_data.size() / size;
      for (int i = 0; i < size; i++) {
        auto file_path =
            dir / (std::string(algname) + "_" + std::to_string(pop) + "_" +
                   std::string(proname) + "_F" + std::to_string(i + 1) + "_" +
                   std::to_string(dim) + "_curve.txt");
        std::ofstream os(file_path.string());
        for (int j = 0; j < len; j++)
          os << graph_data[i * len + j] << '\n';
      }
    }
  }
  template <bool b> void reserve(int sz) {
    if constexpr (b) {
      table_data.reserve(sz);
    } else {
      graph_data.reserve(sz);
    }
  }

  template <bool b, typename V> auto insert(const V &v) {
    if constexpr (b)
      table_data.insert(table_data.end(), std::begin(v), std::end(v));
    else
      graph_data.insert(graph_data.end(), std::begin(v), std::end(v));
  }
};

template <typename V>
concept benchmark_result_concept = requires(V &v) {
  typename V::value_type;
  v.template insert<false>(std::vector<typename V::value_type>());
  v.template insert<true>(std::vector<typename V::value_type>());
};

template <std::uint64_t Alg_HashName, std::uint64_t Prob_HashName,
          bool Memory_Flag, int Dim, int Pop_Size, int Max, int Runs, int P = 1,
          typename V,
          typename Parameter_Type = typename single_algorithm_parameter<
              Alg_HashName, typename V::value_type>::base_type>
  requires benchmark_result_concept<V> &&
           single_algorithm_problem_concept<
               Alg_HashName, Prob_HashName, Dim, Pop_Size, Max, Memory_Flag, 1,
               Parameter_Type, typename V::value_type> &&
           (Runs >= 1)
inline auto single_benchmark(
    V &v, const Parameter_Type &par = Parameter_Type()) noexcept {
  using T = typename V::value_type;
  using F = single_problem<Prob_HashName, P, Dim, T>;
  std::array<T, Runs> result;
  std::conditional_t<Memory_Flag, tool::curve_t<true, T, Max>, void *> data{};
  auto t0 = std::chrono::high_resolution_clock::now();
  using Algorithm =
      single_algorithm<Alg_HashName, Dim, Pop_Size, Max, Memory_Flag>;
  for (int i = 0; i < Runs; i++) {
    if constexpr (!Memory_Flag) {
      auto [_, u] = Algorithm()(F(), F().L, F().U, par);
      result[i] = u;
    } else {
      constexpr T k = T(1) / T(Runs);
      auto [_, u, cur] = Algorithm()(F(), F().L, F().U, par);
      result[i] = u;
      for (int j = 0; j < 2 * Max; j++)
        data[j] += cur[j] * k;
    }
  }
  if constexpr (Memory_Flag) {
    v.template insert<false>(data);
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  auto [m, s, b, b1, b2, b3, b4] =
      tool::mean_std((T *)(result.data()), (T *)(result.data() + Runs));
  std::array<T, 8> box_table;
  box_table[0] = m;
  box_table[1] = s;
  box_table[2] = b;
  box_table[3] =
      static_cast<T>(
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
              .count()) /
      T(Runs);
  box_table[4] = b1;
  box_table[5] = b2;
  box_table[6] = b3;
  box_table[7] = b4;
  v.template insert<true>(box_table);
  if constexpr (P < F::size)
    single_benchmark<Alg_HashName, Prob_HashName, Memory_Flag, Dim, Pop_Size,
                     Max, Runs, P + 1>(v, par);
}

template <std::floating_point T, std::uint64_t Alg_HashName,
          std::uint64_t Prob_HashName, bool Memory_Flag, int Dim, int Pop_Size,
          int Max, int Runs, int P = 1, typename V, typename Task,
          typename Parameter_Type =
              typename single_algorithm_parameter<Alg_HashName, T>::base_type>
  requires single_algorithm_problem_task_concept<
               Alg_HashName, Prob_HashName, Dim, Pop_Size, Max, Memory_Flag, P,
               Parameter_Type, T, Task> &&
           (Runs >= 1)
inline auto parallel_benchmark(
    V &v, Task &pool, const Parameter_Type &par = Parameter_Type()) noexcept {
  using F = single_problem<Prob_HashName, P, Dim, T>;
  using Algorithm =
      single_algorithm<Alg_HashName, Dim, Pop_Size, Max, Memory_Flag>;
  for (int i = 0; i < Runs; i++)
    v[(P - 1) * Runs + i] = pool.submit(
        benchamrk_calc<F, Algorithm, Memory_Flag, Parameter_Type>, par);
  if constexpr (P < F::size)
    parallel_benchmark<T, Alg_HashName, Prob_HashName, Memory_Flag, Dim,
                       Pop_Size, Max, Runs, P + 1>(v, pool, par);
}
template <std::uint64_t Alg_HashName, std::uint64_t Prob_HashName,
          bool Memory_Flag, int Dim, int Pop_Size, int Max, int Runs,
          bool parallel, typename V, typename Task,
          typename Parameter_Type = typename single_algorithm_parameter<
              Alg_HashName, typename V::value_type>::base_type>
  requires benchmark_result_concept<V> &&
           ((parallel &&
             single_algorithm_problem_task_concept<
                 Alg_HashName, Prob_HashName, Dim, Pop_Size, Max, Memory_Flag,
                 1, Parameter_Type, typename V::value_type, Task>) ||
            (!parallel &&
             single_algorithm_problem_concept<
                 Alg_HashName, Prob_HashName, Dim, Pop_Size, Max, Memory_Flag,
                 1, Parameter_Type, typename V::value_type>)) &&
           (Runs >= 1)
inline auto template_benchamrk(
    V &result, Task &pool,
    const Parameter_Type &par = Parameter_Type()) noexcept {
  using T = typename V::value_type;
  using F = single_problem<Prob_HashName, 1, Dim, T>;
  using Algorithm =
      single_algorithm<Alg_HashName, Dim, Pop_Size, Max, Memory_Flag>;
  if constexpr (!parallel) {
    single_benchmark<Alg_HashName, Prob_HashName, Memory_Flag, Dim, Pop_Size,
                     Max, Runs>(result, par);
  } else {
    using Type = std::conditional_t<
        Memory_Flag,
        std::pair<T, typename std::tuple_element_t<Memory_Flag ? 2 : 0,
                                                   decltype(Algorithm()(
                                                       F(), F().L, F().U))>>,
        T>;
    std::vector<std::future<Type>> f(F::size * Runs);
    parallel_benchmark<T, Alg_HashName, Prob_HashName, Memory_Flag, Dim,
                       Pop_Size, Max, Runs>(f, pool, par);
    for (int i = 0; i < F::size; i++) {
      std::array<T, Runs> tmp_result{};
      std::conditional_t<Memory_Flag, std::array<T, 2 * Max>, void *>
          tmp_data{};
      auto t0 = std::chrono::high_resolution_clock::now();
      for (int j = 0; j < Runs; j++) {
        if constexpr (!Memory_Flag) {
          tmp_result[j] = f[i * Runs + j].get();
        } else {
          auto [u, v] = f[i * Runs + j].get();
          tmp_result[j] = u;
          for (int k = 0; k < 2 * Max; k++)
            tmp_data[k] += v[k] / T(Runs);
        }
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      if constexpr (Memory_Flag) {
        result.template insert<false>(tmp_data);
      }
      auto [m, s, b, b1, b2, b3, b4] = tool::mean_std(
          (T *)(tmp_result.data()), (T *)(tmp_result.data() + Runs));
      std::array<T, 8> box_table;
      box_table[0] = m;
      box_table[1] = s;
      box_table[2] = b;
      box_table[3] =
          static_cast<T>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                  .count()) /
          T(Runs);
      box_table[4] = b1;
      box_table[5] = b2;
      box_table[6] = b3;
      box_table[7] = b4;
      result.template insert<true>(box_table);
    }
  }
}
} // namespace sevobench_detail

template <std::uint64_t Alg_HashName, std::uint64_t Prob_HashName,
          bool Memory_Flag = false, int Dim = 30, int Pop_Size = 30,
          int Max = (Memory_Flag ? 1000 * Dim / Pop_Size : 1000 * Dim),
          int Runs = 30, bool parallel = true, std::floating_point T = float,
          typename Parameter_Type =
              typename single_algorithm_parameter<Alg_HashName, T>::base_type>
  requires sevobench_detail::single_algorithm_problem_task_concept<
               Alg_HashName, Prob_HashName, Dim, Pop_Size, Max, Memory_Flag, 1,
               Parameter_Type, T, parallel_task> &&
           (Runs >= 1)
inline auto evo_bench(const Parameter_Type &par = Parameter_Type()) noexcept {
  using Algorithm =
      single_algorithm<Alg_HashName, Dim, Pop_Size, Max, Memory_Flag>;
  using F = single_problem<Prob_HashName, 1, Dim, T>;
  sevobench_detail::benchmark_result<T, Algorithm, F> result{};
  result.template reserve<true>(8 * F::size);
  if (!init_single_problem<Prob_HashName, Dim, T>::is_init) {
    init_single_problem<Prob_HashName, Dim, T>{};
  }
  if constexpr (Memory_Flag)
    result.template reserve<false>(2 * Max * F::size);
  if constexpr (parallel) {
    parallel_task pt{std::min(std::thread::hardware_concurrency(),
                              (unsigned int)(Runs * F::size))};
    sevobench_detail::template_benchamrk<Alg_HashName, Prob_HashName,
                                         Memory_Flag, Dim, Pop_Size, Max, Runs,
                                         parallel>(result, pt, par);
  } else {
    sevobench_detail::single_benchmark<Alg_HashName, Prob_HashName, Memory_Flag,
                                       Dim, Pop_Size, Max, Runs>(result, par);
  }
  return result;
}

} // namespace sevobench
