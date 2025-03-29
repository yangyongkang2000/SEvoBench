#include "SEvoBench/sevobench.hpp"
#include <iostream>
template <int N, typename T> T sphere(const T *x) {
  return std::inner_product(x, x + N, x, T(0));
}
template <std::uint64_t Alg_Hash, std::floating_point T, int Dim = 30,
          int Pop_Size = 100, int Max = 10000 * Dim>
void test_algorithm() {
  auto [u, v] =
      sevobench::other_algorithm::single_algorithm<Alg_Hash, Dim, Pop_Size, Max,
                                                   false>()(sphere<Dim, T>,
                                                            T(-100), T(100));
  auto [u1, v1, w1] =
      sevobench::other_algorithm::single_algorithm<Alg_Hash, Dim, Pop_Size,
                                                   Max / Pop_Size, true>()(
          sphere<Dim, T>, T(-100), T(100));
  if (!(v >= 0 && v < T(1e-3))) {
    std::cout << "ALGORITHM IS FAILED!\n";
  }
  if (!(v1 >= 0 && v1 < T(1e-3))) {
    std::cout << "ALGORITHM IS FAILED!\n";
  }
  if (!(std::abs(sphere<Dim>(u.data()) - v) < T(1e-6))) {
    std::cout << "ALGORITHM IS FAILED!\n";
  }
  if (!(std::abs(sphere<Dim>(u1.data()) - v1) < T(1e-6))) {
    std::cout << "ALGORITHM IS FAILED!\n";
  }
  for (int i = 0; i < (Max / Pop_Size - 1); i++) {
    if (!(w1[2 * i] < w1[2 * (i + 1)])) {
      std::cout << "ALGORITHM IS FAILED!\n";
    }
    if (!(w1[2 * i + 1] >= w1[2 * (i + 1) + 1])) {
      std::cout << "ALGORITHM IS FAILED!\n";
    }
  }
}
template <int Pop_Size = 100, int Max = 300> void test_lshade() {
  if (!((sevobench::other_algorithm::shade_detail::lshade_ite<Pop_Size, 4>(
            sevobench::other_algorithm::shade_detail::iteration_maxfes<
                Pop_Size, 4, Max>)) == Max)) {
    std::cout << "ALGORITHM IS FAILED!\n";
  }
}

#define TEST_ALGORITHM(ALG)                                                    \
  test_algorithm<ALG, float>();                                                \
  test_algorithm<ALG, double>();                                               \
  test_algorithm<ALG, long double>();
int main() {
  using namespace sevobench::other_algorithm;
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
