
#include "SEvoBench/sevobench.hpp"
#include <cassert>
#include <iostream>
#include <set>

template <int N> void test_parallel_task() {
  auto sz = std::thread::hardware_concurrency();
  sevobench::parallel_task pool(sz);
  std::vector<int> ids(N * sz);
  std::vector<std::future<void>> f(N * sz);
  for (decltype(sz) i = 0; i < N * sz; i++) {
    f[i] = pool.submit([&](auto i) { ids[i] = i; }, i);
  }
  for (decltype(sz) i = 0; i < N * sz; i++) {
    f[i].get();
    assert(ids[i] == i);
  }
}

int main() { test_parallel_task<1>(); }