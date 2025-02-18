#include "SEvoBench/sevobench.hpp"
template <int N, typename T> T sphere(const T *x) {
  return std::inner_product(x, x + N, x, T(0));
}
int main() {
  constexpr int Dim = 10;
  constexpr int Pop_Size = 30;
  using T = float;
  std::vector<std::vector<T>> pos(Pop_Size);
  std::begin(pos);
  for (auto &_ : pos)
    _.resize(Dim);
  sevobench::other_algorithm::de_optimize<Dim, Pop_Size>(pos, sphere<Dim, T>,
                                                         T(-1), T(1));
  sevobench::other_algorithm::abc_optimize<Dim, Pop_Size>(pos, sphere<Dim, T>,
                                                          T(-1), T(1));
  sevobench::other_algorithm::pso_optimize<Dim, Pop_Size>(pos, sphere<Dim, T>,
                                                          T(-1), T(1));
  sevobench::other_algorithm::es_optimize<Dim, Pop_Size>(pos[0], sphere<Dim, T>,
                                                         T(-1), T(1));
  sevobench::other_algorithm::random_search_optimize<Dim, Pop_Size>(
      pos[0], sphere<Dim, T>, T(-1), T(1));
  sevobench::other_algorithm::template_pso_optimize<true, Dim, Pop_Size, 1000,
                                                    false>(pos, sphere<Dim, T>,
                                                           T(-1), T(1));
  sevobench::other_algorithm::template_pso_optimize<false, Dim, Pop_Size, 1000,
                                                    false>(pos, sphere<Dim, T>,
                                                           T(-1), T(1));
  sevobench::other_algorithm::slpso_optimize<Dim, Pop_Size, 1000>(
      pos, sphere<Dim, T>, T(-1), T(1));
  sevobench::other_algorithm::cso_optimize<Dim, Pop_Size>(pos, sphere<Dim, T>,
                                                          T(-1), T(1));
  sevobench::other_algorithm::jade_optimize<Dim, Pop_Size>(pos, sphere<Dim, T>,
                                                           T(-1), T(1));
  sevobench::other_algorithm::template_shade_optimize<true, Dim, Pop_Size, 1000,
                                                      false>(
      pos, sphere<Dim, T>, T(-1), T(1));
  sevobench::other_algorithm::template_shade_optimize<false, Dim, Pop_Size,
                                                      1000, false>(
      pos, sphere<Dim, T>, T(-1), T(1));
}