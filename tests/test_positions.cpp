//
// Created by 杨永康 on 2024/2/24.
//
#include"SEvoBench/single_algorithm.hpp"
#include<valarray>
template <int N, typename T> T sphere(const T *x) {
  return std::inner_product(x, x + N, x, T(0));
}
int main() {
  constexpr int Dim=10;
  constexpr int Pop_Size=30;
  using T=float;
  std::vector<std::vector<T>> pos(Pop_Size);
  std::begin(pos);
  for(auto &_:pos)
    _.resize(Dim);
  sevobench::de_optimize<Dim,Pop_Size>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::abc_optimize<Dim,Pop_Size>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::pso_optimize<Dim,Pop_Size>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::es_optimize<Dim,Pop_Size>(pos[0],sphere<Dim,T>,T(-1),T(1));
  sevobench::random_search_optimize<Dim,Pop_Size>(pos[0],sphere<Dim,T>,T(-1),T(1));
  sevobench::template_pso_optimize<true,Dim,Pop_Size,1000,false>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::template_pso_optimize<false,Dim,Pop_Size,1000,false>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::slpso_optimize<Dim,Pop_Size,1000>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::cso_optimize<Dim,Pop_Size>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::jade_optimize<Dim,Pop_Size>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::template_shade_optimize<true,Dim,Pop_Size,1000,false>(pos,sphere<Dim,T>,T(-1),T(1));
  sevobench::template_shade_optimize<false,Dim,Pop_Size,1000,false>(pos,sphere<Dim,T>,T(-1),T(1));
}