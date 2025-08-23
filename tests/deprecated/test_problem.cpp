
#include "SEvoBench/single_problem.hpp"
#include <bitset>
#include <cassert>
#include <iostream>
template <typename T, int N, typename M>
inline void test_orthogonal(const M &mat) noexcept {
  std::vector<std::array<T, N>> product(N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      product[i][j] = std::inner_product(mat[i].begin(), mat[i].end(),
                                         mat[j].begin(), T(0));
      assert(std::abs(product[i][j] - T(i == j ? 1 : 0)) < T(1e-1));
    }
  }
}

template <std::uint64_t P, int Dim, int I1, int I2, typename T, typename V,
          typename G, typename H>
auto test_problem_data(const V &v, const G &g, const H &h) {

  using Prob = sevobench::single_problem<P, I1, Dim, T>;
  static_assert(Prob::index == I1 && Prob::dim == Dim);
  assert(std::all_of(g[I1 - 1].begin(), g[I1 - 1].end(),
                     [](auto &x) { return x >= Prob::L && x <= Prob::U; }));
  test_orthogonal<T, Dim>(h[I1 - 1]);
  assert(std::abs(Prob()(g[I1 - 1].data()) - v[I1 - 1]) < T(1e-1));
  if constexpr (I1 < I2)
    return test_problem_data<P, Dim, I1 + 1, I2, T>(v, g, h);
}

template <std::uint64_t P, int Dim, int S, int I1, int I2, typename T,
          typename V, typename G, typename H, typename... Arg>
auto test_extra_problem_data(const V &v, const G &g, const H &h,
                             const Arg &...arg) {
  using Prob = sevobench::single_problem<P, I1, Dim, T>;
  static_assert(Prob::index == I1 && Prob::dim == Dim);
  auto len = g.size() / Dim;
  for (decltype(len) i = 0; i < len; i++) {
    assert(std::all_of(g.begin() + i * Dim, g.begin() + (i + 1) * Dim,
                       [](auto x) { return x >= Prob::L && x <= Prob::U; }));
    assert((std::abs(Prob()(g.data() + i * Dim) - v[I1 - S][i]) < T(1e-1)));
  }
  for (auto &_ : h)
    test_orthogonal<T, Dim>(_);
  if constexpr (I1 < I2)
    return test_extra_problem_data<P, Dim, S, I1 + 1, I2, T>(v, arg...);
}

template <int Dim, typename V> void test_perm_data(const V &v) {
  for (auto &tmp : v) {
    std::bitset<Dim> d;
    for (auto &_ : tmp)
      d.set(_);
    assert(d.all());
  }
}

template <int P, int Dim, typename T> void test_problem() {

  if constexpr (P == 0) {
    assert(((sevobench::yao_func::detail::is_init_rotate<T, Dim>) == false) &&
           ((sevobench::yao_func::detail::is_init_shift<T, Dim>) == false));
    sevobench::init_single_problem<sevobench::ShiftFunc, Dim, T>{};
    assert(((sevobench::yao_func::detail::is_init_rotate<T, Dim>) == false) &&
           ((sevobench::yao_func::detail::is_init_shift<T, Dim>) == true));
    sevobench::init_single_problem<sevobench::AsyRotateShiftFunc, Dim, T>{};
    assert(((sevobench::yao_func::detail::is_init_rotate<T, Dim>) == true) &&
           ((sevobench::yao_func::detail::is_init_shift<T, Dim>) == true));
    sevobench::init_single_problem<sevobench::RotateShiftFunc, Dim, T>{};
    sevobench::init_single_problem<sevobench::AsyShiftFunc, Dim, T>{};
    constexpr T v[13] = {};
    test_problem_data<sevobench::ShiftFunc, Dim, 1, 13, T>(
        v, sevobench::yao_func::detail::shift_vector_data<T, Dim>,
        sevobench::yao_func::detail::rotate_matrix_data<T, Dim>);
    test_problem_data<sevobench::RotateShiftFunc, Dim, 1, 13, T>(
        v, sevobench::yao_func::detail::shift_vector_data<T, Dim>,
        sevobench::yao_func::detail::rotate_matrix_data<T, Dim>);
    test_problem_data<sevobench::AsyShiftFunc, Dim, 1, 13, T>(
        v, sevobench::yao_func::detail::shift_vector_data<T, Dim>,
        sevobench::yao_func::detail::rotate_matrix_data<T, Dim>);
    test_problem_data<sevobench::AsyRotateShiftFunc, Dim, 1, 13, T>(
        v, sevobench::yao_func::detail::shift_vector_data<T, Dim>,
        sevobench::yao_func::detail::rotate_matrix_data<T, Dim>);
  } else {
    constexpr auto Pr = P == 1 ? sevobench::CEC2020 : sevobench::CEC2022;
    constexpr T v[] = {100,  1100, 700,  1900, 1700,
                       1600, 2100, 2200, 2400, 2500};
    constexpr T v1[] = {300,  400,  600,  800,  900,  1800,
                        2000, 2200, 2300, 2400, 2600, 2700};
    using Init = sevobench::init_single_problem<Pr, Dim, T>;
    assert(Init::is_init == false);
    Init{};
    assert(Init::is_init == true);
    if constexpr (P == 1) {
      test_problem_data<Pr, Dim, 1, 7, T>(
          v, sevobench::ieee_cec_set::cec2020_data::shift_vector_data<T, Dim>,
          sevobench::ieee_cec_set::cec2020_data::rotate_matrix_data<T, Dim>);
      const std::vector<std::vector<T>> extra{{2200, 2300, 2400},
                                              {2400, 2500, 2600, 2700},
                                              {2500, 2600, 2700, 2800, 2900}};
      test_extra_problem_data<Pr, Dim, 8, 8, 10, T>(
          extra, sevobench::ieee_cec_set::cec2020_data::shift_o8<T, Dim>,
          sevobench::ieee_cec_set::cec2020_data::mat8<T, Dim>,
          sevobench::ieee_cec_set::cec2020_data::shift_o9<T, Dim>,
          sevobench::ieee_cec_set::cec2020_data::mat9<T, Dim>,
          sevobench::ieee_cec_set::cec2020_data::shift_o10<T, Dim>,
          sevobench::ieee_cec_set::cec2020_data::mat10<T, Dim>);
      test_perm_data<Dim>(
          sevobench::ieee_cec_set::cec2020_data::perm_data<Dim>);
    } else {
      test_problem_data<Pr, Dim, 1, 8, T>(
          v1, sevobench::ieee_cec_set::cec2022_data::shift_vector_data<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::rotate_matrix_data<T, Dim>);
      const std::vector<std::vector<T>> extra{
          {2300, 2500, 2600, 2400, 2700},
          {2400, 2600, 2500},
          {2600, 2800, 2900, 3000, 2800},
          {2700, 3000, 3200, 2800, 3100, 2900}};
      test_extra_problem_data<Pr, Dim, 9, 9, 12, T>(
          extra, sevobench::ieee_cec_set::cec2022_data::shift_o9<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::mat9<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::shift_o10<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::mat10<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::shift_o11<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::mat11<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::shift_o12<T, Dim>,
          sevobench::ieee_cec_set::cec2022_data::mat12<T, Dim>);
      test_perm_data<Dim>(
          sevobench::ieee_cec_set::cec2022_data::perm_data<Dim>);
    }
  }
}

#define TEST_PROBLEMS(Dim, T)                                                  \
  test_problem<0, Dim, T>();                                                   \
  test_problem<1, Dim, T>();                                                   \
  test_problem<2, Dim, T>();
#define TEST_PROBLEM(Dim)                                                      \
  TEST_PROBLEMS(Dim, float)                                                    \
  TEST_PROBLEMS(Dim, double) TEST_PROBLEMS(Dim, long double)

int main() {
  TEST_PROBLEM(1)
  TEST_PROBLEM(2)
  TEST_PROBLEM(4)
  TEST_PROBLEM(8)
  TEST_PROBLEM(10)
  TEST_PROBLEM(20)
  TEST_PROBLEM(30)
  TEST_PROBLEM(40)
  TEST_PROBLEM(50)
  TEST_PROBLEM(60)
  TEST_PROBLEM(70)
  TEST_PROBLEM(80)
  TEST_PROBLEM(90)
  TEST_PROBLEM(100)
}