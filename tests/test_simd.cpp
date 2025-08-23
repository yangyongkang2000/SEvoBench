#include "SEvoBench/problem/large_scale/cec2010.hpp"
#include "SEvoBench/sevobench.hpp"
#include <iostream>
template <int Dim, std::floating_point T> void test() {
  using namespace sevobench;
  std::array<T, Dim> x;
  std::random_device rd{};
  std::default_random_engine gen{rd()};
  std::uniform_real_distribution<T> dis(T(-10), T(10));
  std::generate_n(x.data(), Dim, [&] { return dis(gen); });
  {
    auto v1 =
        problem::simd::ackley<true, Dim, problem::simd_width<T>, T>(x.data());
    auto v2 =
        problem::simd::ackley<false, Dim, problem::simd_width<T>, T>(x.data());
    std::cout << v1 << "," << v2 << '\n';
  }
  {
    auto v1 = problem::simd::rastrigin<true, Dim, problem::simd_width<T>, T>(
        x.data());
    auto v2 = problem::simd::rastrigin<false, Dim, problem::simd_width<T>, T>(
        x.data());
    std::cout << v1 << "," << v2 << '\n';
  }
  {
    auto v1 =
        problem::simd::elliptic<true, Dim, problem::simd_width<T>, T>(x.data());
    auto v2 = problem::simd::elliptic<false, Dim, problem::simd_width<T>, T>(
        x.data());
    std::cout << v1 << "," << v2 << '\n';
  }
}
int main() { test<39, float>(); }