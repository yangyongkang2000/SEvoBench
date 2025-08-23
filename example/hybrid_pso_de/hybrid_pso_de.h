#pragma once
#include "SEvoBench/sevobench.hpp"
template <std::floating_point T>
void hybrid_pso_de(auto &&pso_alg, auto &&de_alg, auto &&pop, auto &&f, T lb,
                   T ub, auto &&alg) {
  using namespace sevobench;
  auto ps = pop.pop_size();
  auto dim = pop.dim();
  population<T> pop1(ps, dim), pop2(ps, dim), pop3(ps, dim);
  pso_module::pso_velocity<T> vec(ps, std::vector<T>(dim));
  for (auto &_ : pop)
    _.evaluate(f);
  alg.add_fes(pop.pop_size());
  pso_alg.topology()->prepare(pop);
  do {
    pop1 = pop;
    pop3 = pop;
    pso_alg.iterator(pop1, vec, f, lb, ub, alg);
    de_alg.iterator(pop3, pop2, f, lb, ub, alg);
    for (int i = 0; i < ps; i++)
      for (int j = 0; j < dim; j++)
        vec[i][j] = pop3[i][j] - pop[i][j];
    for (int i = 0; i < ps; i++)
      std::swap(pop[i],
                pop1[i].fitness() < pop3[i].fitness() ? pop1[i] : pop3[i]);
  } while (alg.current_fes() < alg.max_fes());
}
