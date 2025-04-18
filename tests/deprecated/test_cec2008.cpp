#include "SEvoBench/cec2008.hpp"
#include <cassert>
#include <iostream>
template <typename T> std::array<std::array<T, 6>, 2> tmp_result;

template <std::uint64_t P, int Dim, typename T, int Index, int I = 1>
void test_cec2008() {
  using Pro = sevobench::single_problem<P, I, Dim, T>;
  std::array<T, Dim> tmp{};
  assert(
      std::abs(Pro()(sevobench::ieee_cec_set::cec2008_data::shift_vector_data<
                         T, Dim>[I - 1]
                         .data()) -
               (I == 3 ? (Dim - 1) : 0)) < T(1e-3));
  tmp_result<T>[Index][I - 1] = Pro()(tmp.data());
  if constexpr (I < Pro::size)
    test_cec2008<P, Dim, T, Index, I + 1>();
};

template <std::uint64_t P, int Dim, typename T, int Index> struct test_op {
  test_op() {
    std::cout << sevobench::simd_id() << "\n";
    test_cec2008<P, Dim, T, Index>();
  }
};

auto init1 = sevobench::init_single_problem<sevobench::CEC2008, 512, float>{};
auto init2 = sevobench::init_single_problem<sevobench::CEC2008, 512, double>{};

#define TEST(P, Dim, T, Index, p) auto p = test_op<P, Dim, T, Index>();

TEST(sevobench::CEC2008, 512, float, 0, p1)
TEST(sevobench::CEC2008, 512, double, 0, p2)

#undef INSTRSET

TEST(sevobench::CEC2008, 512, float, 1, p3)
TEST(sevobench::CEC2008, 512, double, 1, p4)

int main() {
  for (int i = 0; i < 6; i++) {
    assert(std::abs(tmp_result<double>[0][i] - tmp_result<double>[1][i]) <
           double(1e-3));
    assert(std::abs(tmp_result<float>[0][i] - tmp_result<float>[1][i]) <
           double(1e-3));
  }
}