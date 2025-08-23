#include "SEvoBench/sevobench.hpp"
#include <chrono>
#include <iostream>
template<int N, typename T>
inline auto real_func(std::span<const T> t) noexcept {
    T sum = 0;
    for (int i = 0; i < N - 1; i++)
        sum += sevobench::tool::Pow<2>(t[i] - 1) +
               100 * sevobench::tool::Pow<2>((t[i + 1] - t[i] * t[i]));
    return sum;
}
template<std::floating_point T>
auto test_shade() {
    using namespace sevobench::de_module;
    constexpr int Pop_Size = 100;
    constexpr int Dim = 30;
    sevobench::evolutionary_algorithm alg(Dim * 1000);
    sevobench::population<T> pop(Pop_Size, Dim, T(-100), T(100));
    auto de =
            de_algorithm_builder<T>()
                    .mutation(
                            std::make_unique<ttpb1_mutation<T>>())
                    .parameter(std::make_unique<shade_parameter<T>>())
                    .constraint_handler(std::make_unique<midpoint_target_repair<T>>())
                    .crossover(std::make_unique<binomial_crossover<T>>())
                    .archive(std::make_unique<random_archive<T>>());
    de.population_strategy(std::make_unique<de_population<T>>())
            .build()
            .run(pop, real_func<Dim, T>, T(-100), T(100), alg);
    return std::min_element(
                   pop.begin(), pop.end(),
                   [](auto &x, auto &y) { return x.fitness() < y.fitness(); })
            ->fitness();
}
template<std::floating_point T>
std::tuple<T, T, double> analyze_runs(int runs) {
    std::vector<T> results;
    std::vector<double> timings;
    results.reserve(runs);
    timings.reserve(runs);

    for (int i = 0; i < runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        T value = test_shade<T>();
        auto end = std::chrono::high_resolution_clock::now();

        results.push_back(value);
        timings.push_back(
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                        .count());
    }

    auto result =
            sevobench::tool::mean_std(results.data(), results.data() + runs);
    double avg_time = std::accumulate(timings.begin(), timings.end(), 0.0) / runs;

    return {result[0], result[1], avg_time};
}
int main() {
    using T = double;
    constexpr int Runs = 10;
    analyze_runs<T>(Runs);
    auto [mean, variance, avg_ms] = analyze_runs<double>(Runs);
    std::cout << "Results after " << Runs << " runs:\n"
              << "  Average value: " << mean << "\n"
              << "  Variance:      " << variance << "\n"
              << "  Average time:   " << avg_ms << " ms\n";
    return 0;
}
