#include <iostream>
#include"hybrid_pso_de.h"
template <typename T>
inline auto rosenbrock(std::span<const T> t) noexcept {
    T sum = 0;
    for (int i = 0; i < t.size() - 1; i++)
        sum += sevobench::tool::Pow<2>(t[i] - 1) +
               100 * sevobench::tool::Pow<2>((t[i + 1] - t[i] * t[i]));
    return sum;
}
int main() {
    using namespace sevobench;
    using namespace sevobench::pso_module;
    using namespace sevobench::de_module;
    auto topo=std::make_unique<lbest_topology<float>>();
    auto upda= std::make_unique<spherical_update<float>>();
    auto pso=pso_algorithm_builder().
            topology(std::move(topo)).
            update(std::move(upda)).
            build();
    auto p= std::make_unique<shade_parameter<float>>();
    auto m=std::make_unique<ttpb1_mutation<float>>();
    auto c= std::make_unique<binomial_crossover<float>>();
    auto h=std::make_unique<midpoint_target_repair<float>>();
    auto s= std::make_unique<de_population<float>>();
    auto shade=de_algorithm_builder().
            mutation(std::move(m)).
            parameter(std::move(p)).
            population_strategy(std::move(s)).
            crossover(std::move(c)).
            constraint_handler(std::move(h)).
            build();
    constexpr int dim=30;
    constexpr int pop_size=100;
    constexpr int max_fes=1000*dim;
    evolutionary_algorithm alg(max_fes,pop_size,dim);
    population<float> pop(pop_size,dim,float(-100),float(100));
    hybrid_pso_de(pso,shade,pop,[](std::span<const float> x){return rosenbrock(x);},float(-100),float(100),alg);
    std::cout << "best value:"<<std::min_element(pop.begin(),pop.end(),[](auto &l,auto &r){
        return l.fitness()<r.fitness();
    })->fitness() << std::endl;
    return 0;
}
