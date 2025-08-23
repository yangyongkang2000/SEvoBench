#include "SEvoBench/sevobench.hpp"
#include <cassert>
#include <iostream>

template <typename T> void test_find_median() {
  using namespace sevobench::tool;
  constexpr T num[6] = {1, 2, 3, 4, 5, 6};
  static_assert(find_median(num, num + 1) == T(1));
  static_assert(find_median(num, num + 2) == T(1.5));
  static_assert(find_median(num, num + 5) == T(3));
  static_assert(find_median(num, num + 6) == T(3.5));
}

template <typename T> void test_mean_std() {
  std::vector<T> data{5, 4, 3, 2, 1};
  auto [mean, std1, min, median_first_half, median, median_second_half, max] =
      sevobench::tool::mean_std(data.begin(), data.end());
  if (mean == 3) {
    std::cout << "mean is true!\n";
  }
  if (std::abs(std1 - std::sqrt(2.5)) < T(1e-6)) {
    std::cout << "std is true!\n";
  }
  if (min == T(1)) {
    std::cout << "min is true!\n";
  }
  if (median_first_half == T(1.5)) {
    std::cout << "median_first_half is true!\n";
  }
  if (median == T(3)) {
    std::cout << "median is true!\n";
  }
  if (median_second_half == T(4.5)) {
    std::cout << "median_second_half is true!\n";
  }
  if (max == T(5)) {
    std::cout << "max is true!\n";
  }
}

void test_kunth_shffule() {
  int num[] = {0, 1, 2, 3, 4};
  std::random_device rd{};
  for (int i = 0; i < 5; i++)
    sevobench::tool::kunth_shuffle(num + i, num + i + 1, rd());
  for (int i = 0; i < 5; i++)
    if (!(num[i] == i)) {
      std::cout << "kunth_shffule is failed!\n";
    }
  sevobench::tool::kunth_shuffle(num, num + 5, rd());
  std::sort(num, num + 5);
  for (int i = 0; i < 5; i++)
    if (!(num[i] == i)) {
      std::cout << "kunth_shffule is failed!\n";
    }
}

#define TEST_STATS(T)                                                          \
  test_mean_std<T>();                                                          \
  test_find_median<T>();
#define TEST_STAT() TEST_STATS(float) TEST_STATS(double) TEST_STATS(long double)

int main() {
  test_kunth_shffule();
  TEST_STAT()
  return 0;
}
