#pragma once
#include "../../common/tool.hpp"
#include "../problem.hpp"
#include "cec_base_problem.hpp"
#include <fstream>
#include <sstream>
namespace sevobench::problem {

namespace cec_detail {

template <std::size_t N, bool only_shift = false, std::floating_point T>
inline void sr_func(std::span<const T> x, std::span<T, N> y,
                    std::span<const T> o, std::span<const T> m) noexcept {
  if constexpr (only_shift) {
    for (std::size_t i = 0; i < N; i++)
      y[i] = x[i] - o[i];
  } else {
    std::array<T, N> tmp;
    for (std::size_t i = 0; i < N; i++)
      tmp[i] = x[i] - o[i];
    for (std::size_t i = 0; i < N; i++)
      y[i] = std::inner_product(tmp.data(), tmp.data() + N, m.data() + i * N,
                                T(0));
  }
}

template <int N, int Num, std::floating_point T>
inline auto cf_cal(std::span<const T> x, std::span<const T> o,
                   std::span<const T, Num> fit,
                   std::span<const T, Num> delta) noexcept {
  std::array<T, Num> w{};
  for (int i = 0; i < Num; i++) {
    for (int j = 0; j < N; j++)
      w[i] += tool::Pow<2>(x[j] - o[i * N + j]);
  }
  auto it = std::find(w.begin(), w.end(), T(0));
  if (it != w.end()) {
    return fit[it - w.begin()];
  }
  for (int i = 0; i < Num; i++)
    w[i] = std::exp(-w[i] / (2 * N * tool::Pow<2>(delta[i]))) / std::sqrt(w[i]);
  auto w_sum = std::accumulate(w.begin(), w.end(), T(0));
  return std::inner_product(w.begin(), w.end(), fit.begin(), T(0)) / w_sum;
}

template <int Dim, typename T>
inline auto read_cec_data(const std::string &file_name) noexcept {
  std::vector<T> result(Dim);
  std::ifstream is(file_name);
  for (int i = 0; i < Dim; i++) {
    is >> result[i];
  }
  return result;
}

template <int Dim, std::floating_point T, std::size_t M, typename F>
  requires(Dim > 0) && (Dim % 10 == 0)
inline auto calculate_hybrid(std::span<T, Dim> y,
                             const std::array<std::pair<int, int>, M> &ratios,
                             const std::array<F, M> &functions) noexcept {
  T sum(0);
  int offset = 0;
  for (std::size_t i = 0; i < M; i++) {
    auto [u, v] = ratios[i];
    auto dim = (u * Dim) / v;
    sum += functions[i](y.subspan(offset, dim));
    offset += dim;
  }
  return sum;
}

template <int Dim, bool b = false, std::floating_point T, std::size_t M,
          typename G>
inline auto calculate_composition(
    std::span<const T> x, const std::array<T, M> &deltas,
    const std::array<T, M> &lams, const std::array<T, M> &biases,
    const std::array<G, M> &functions, std::span<const T> shift,
    std::span<const T> rotate, std::span<const int> shuffle) noexcept {
  std::array<T, M> fits;
  std::array<T, Dim> y;
  for (size_t i = 0; i < M; i++) {
    sr_func<Dim, false, T>(x, y, shift.subspan(i * Dim),
                           rotate.subspan(i * Dim * Dim));
    if constexpr (b) {
      std::array<T, Dim> z = y;
      for (int j = 0; j < Dim; j++)
        y[j] = z[shuffle[i * Dim + j]];
    }
    fits[i] = lams[i] * functions[i](y) + biases[i];
  }
  return cf_cal<Dim, M, T>(x, shift, fits, deltas);
}
} // namespace cec_detail

template <int Index, int Dim, std::floating_point T,
          template <int, int, std::floating_point> class Drived>
class cec_common : public problem_common<Index, Dim, T> {
protected:
  std::vector<T> shift;
  std::vector<T> matrix;
  std::vector<int> shuffle;
  int ins = 0;
  void load_rotate_matrix(const std::string &dir_name) {
    std::ostringstream os;
    os << dir_name << '/' << "M_" << Index << "_D" << Dim << ".txt";
    if constexpr (is_composition_problem()) {
      matrix =
          cec_detail::read_cec_data<Dim * Dim * Drived<Index, Dim, T>::cf_num(),
                                    T>(os.str());
    } else {
      matrix = cec_detail::read_cec_data<Dim * Dim, T>(os.str());
    }
  }
  void load_shift_shuffle(const std::string &dir_name) {
    {
      std::ostringstream os;
      os << dir_name << '/' << "shift_data_" << Index << ".txt";
      if constexpr (is_composition_problem()) {
        shift =
            cec_detail::read_cec_data<Dim * Drived<Index, Dim, T>::cf_num(), T>(
                os.str());
      } else {
        shift = cec_detail::read_cec_data<Dim, T>(os.str());
      }
    }
    if constexpr (is_hybrid_problem() || (requires {
                    requires Drived<Index, Dim, T>::is_hybrid_composition();
                  })) {
      constexpr auto num = [] {
        if constexpr (is_composition_problem())
          return Drived<Index, Dim, T>::cf_num();
        else
          return 1;
      }();
      std::ostringstream os;
      os << dir_name << '/' << "shuffle_data_" << Index << "_D" << Dim
         << ".txt";
      shuffle = cec_detail::read_cec_data<num * Dim, int>(os.str());
      for (auto &x : shuffle)
        --x;
    }
  }

public:
  cec_common() = default;
  cec_common(cec_common &&) = default;
  cec_common(const cec_common &) = default;
  cec_common &operator=(cec_common &&) = default;
  cec_common &operator=(const cec_common &) = default;
  cec_common(int _ins) : ins(_ins) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<T> o1(lower_bound(), upper_bound());
    constexpr auto num = [] {
      if constexpr (is_composition_problem())
        return Drived<Index, Dim, T>::cf_num();
      else
        return 1;
    }();
    shift.resize(num * Dim);
    std::generate_n(shift.begin(), num * Dim, [&] { return o1(gen); });
    if constexpr (!is_only_shift()) {
      matrix.resize(Dim * Dim * num);
      for (int i = 0; i < num; i++) {
        tool::generate_rotate_vector<Dim, T>(
            std::span<T>(matrix.data() + i * Dim * Dim, Dim * Dim));
      }
    }
    if constexpr (is_hybrid_problem() || (requires {
                    requires Drived<Index, Dim, T>::is_hybrid_composition();
                  })) {
      shuffle.resize(num * Dim);
      std::default_random_engine gen1(std::random_device{}());
      for (int i = 0; i < num; i++) {
        std::iota(shuffle.data() + i * Dim, shuffle.data() + (i + 1) * Dim,
                  int(0));
        std::shuffle(shuffle.data() + i * Dim, shuffle.data() + (i + 1) * Dim,
                     gen1);
      }
    }
  }
  cec_common(const std::string &dir_name) {
    load_shift_shuffle(dir_name);
    load_rotate_matrix(dir_name);
    ins = 1;
  }
  constexpr static auto lower_bound() { return T(-100); }
  constexpr static auto upper_bound() { return T(100); }
  auto instance() const noexcept { return ins; }
  void set_instance(int _ins) noexcept { ins = _ins; }
  static constexpr auto has_basic_func() {
    return requires(std::span<T, Dim> x) {
      { Drived<Index, Dim, T>::evaluate(x) } -> std::same_as<T>;
    } || requires { requires Drived<Index, Dim, T>::is_basic(); };
  }
  static constexpr auto has_hybrid_func() {
    return requires(std::span<T, Dim> x) {
      { Drived<Index, Dim, T>::hybrid_evaluate(x) } -> std::same_as<T>;
    } || requires { requires Drived<Index, Dim, T>::is_hybrid(); };
  }
  static constexpr auto has_composition_func() {
    return requires(std::span<const T> x, std::span<const int> y) {
      {
        Drived<Index, Dim, T>::composition_evaluate(x, x, x, y)
      } -> std::same_as<T>;
      { Drived<Index, Dim, T>::cf_num() } -> std::same_as<int>;
      requires Drived<Index, Dim, T>::cf_num() > 1;
    } || requires { requires Drived<Index, Dim, T>::is_composition(); };
  }
  static constexpr auto is_basic_problem() {
    return (has_basic_func()) && (!has_hybrid_func()) &&
           (!has_composition_func());
  }
  static constexpr auto is_hybrid_problem() {
    return (!has_basic_func()) && (has_hybrid_func()) &&
           (!has_composition_func());
  }
  static constexpr auto is_composition_problem() {
    return (!has_basic_func()) && (!has_hybrid_func()) &&
           (has_composition_func());
  }
  static constexpr auto is_only_shift() {
    return is_basic_problem() &&
           requires { requires Drived<Index, Dim, T>::is_only_shift(); };
  }
  auto self() const noexcept {
    static_assert(std::is_base_of_v<cec_common<Index, Dim, T, Drived>,
                                    Drived<Index, Dim, T>>,
                  "CRTP FAILED");
    return static_cast<Drived<Index, Dim, T> *>(this);
  }

  static constexpr auto optimum_num() noexcept {
    if constexpr (requires {
                    { Drived<Index, Dim, T>::optimum_num() } -> std::same_as<T>;
                  }) {
      return Drived<Index, Dim, T>::optimum_num();
    } else {
      return 100 * Index;
    }
  }

  auto optimum_solution() const noexcept {
    solution<T> x(shift.data(), shift.data() + Dim);
    x.set_fitness(optimum_num());
    return x;
  }
  auto problem_information() const noexcept {
    return problem_info<T>{.index = Index,
                           .instance = ins,
                           .dim = Dim,
                           .lb = T(-100),
                           .ub = T(100),
                           .optimum = optimum_solution()};
  }
  auto operator()(std::span<const T> x) const noexcept {
    static_assert(is_composition_problem() || is_hybrid_problem() ||
                      is_basic_problem(),
                  "CRTP FAILED");
    if constexpr (!is_composition_problem()) {
      std::array<T, Dim> y;
      cec_detail::sr_func<Dim, is_only_shift(), T>(x, y, shift, matrix);
      if constexpr (is_hybrid_problem()) {
        std::array<T, Dim> z;
        for (int i = 0; i < Dim; i++)
          z[i] = y[shuffle[i]];
        return Drived<Index, Dim, T>::hybrid_evaluate(z) + optimum_num();
      } else {
        return Drived<Index, Dim, T>::evaluate(y) + optimum_num();
      }
    } else {
      return Drived<Index, Dim, T>::composition_evaluate(x, shift, matrix,
                                                         shuffle) +
             optimum_num();
    }
  }
};

} // namespace sevobench::problem
