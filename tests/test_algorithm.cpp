#include "SEvoBench/single_algorithm.hpp"
#include <cassert>
template <int N, typename T> T sphere(const T *x) {
  return std::inner_product(x, x + N, x, T(0));
}
template <std::uint64_t Alg_Hash, std::floating_point T, int Dim = 30,
          int Pop_Size = 100, int Max = 10000 * Dim>
void test_algorithm() {
  auto [u, v] =
      sevobench::single_algorithm<Alg_Hash, Dim, Pop_Size, Max, false>()(
          sphere<Dim, T>, T(-100), T(100));
  auto [u1, v1, w1] =
      sevobench::single_algorithm<Alg_Hash, Dim, Pop_Size, Max / Pop_Size,
                                  true>()(sphere<Dim, T>, T(-100), T(100));
  assert(v >= 0 && v < T(1e-3));
  assert(v1 >= 0 && v1 < T(1e-3));
  assert(std::abs(sphere<Dim>(u.data()) - v) < T(1e-6));
  assert(std::abs(sphere<Dim>(u1.data()) - v1) < T(1e-6));
  for (int i = 0; i < (Max / Pop_Size - 1); i++) {
    assert(w1[2 * i] < w1[2 * (i + 1)]);
    assert(w1[2 * i + 1] >= w1[2 * (i + 1) + 1]);
  }
}
template <int Pop_Size = 100, int Max = 300> void test_lshade() {
  assert(((sevobench::shade_detail::lshade_ite<Pop_Size, 4>(
              sevobench::shade_detail::iteration_maxfes<Pop_Size, 4, Max>)) ==
          Max));
}

#define TEST_ALGORITHM(ALG)                                                    \
  test_algorithm<ALG, float>();                                                \
  test_algorithm<ALG, double>();                                               \
  test_algorithm<ALG, long double>();
int main() {
  using namespace sevobench;
  test_lshade<>();
  TEST_ALGORITHM(PSO)
  TEST_ALGORITHM(SLPSO)
  TEST_ALGORITHM(SPSO2007)
  TEST_ALGORITHM(SPSO2011)
  TEST_ALGORITHM(CSO)
  TEST_ALGORITHM(DE)
  TEST_ALGORITHM(JADE)
  TEST_ALGORITHM(SHADE)
  TEST_ALGORITHM(ABC)
  TEST_ALGORITHM(ES)
}
