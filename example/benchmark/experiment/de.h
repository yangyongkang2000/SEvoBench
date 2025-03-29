#pragma once



#include <vector>
#include <random>
#include <algorithm>
#include <tuple>

template<typename T, typename F>
auto de(F&& f, size_t dim, T lb, T ub, size_t max_evaluations = 10000,
        size_t population_size = 50, T F_factor = 0.8, T CR_prob = 0.9)
-> std::tuple<std::vector<T>, T, size_t>
{
    std::random_device rd;
    std::default_random_engine gen{rd()};
    std::uniform_real_distribution<T> uniform(0.0, 1.0);
    std::uniform_int_distribution<size_t> idx_dist(0, population_size-1);

    std::vector<std::vector<T>> population(population_size, std::vector<T>(dim));
    for (auto& ind : population) {
        std::generate(ind.begin(), ind.end(), [&]{
            return std::uniform_real_distribution<T>(lb, ub)(gen);
        });
    }

    std::vector<std::vector<T>> trials(population_size, std::vector<T>(dim));
    std::vector<T> fitness(population_size), trial_fitness(population_size);

    size_t evaluation_counter = 0;
    for (size_t i = 0; i < population_size; ++i) {
        fitness[i] = f(population[i]);
    }
    evaluation_counter += population_size;

    auto best_it = std::min_element(fitness.begin(), fitness.end());
    auto best_sol = population[std::distance(fitness.begin(), best_it)];
    auto best_fit = *best_it;

    while (evaluation_counter < max_evaluations) {
        for (size_t target = 0; target < population_size; ++target) {
            size_t a, b, c;
            do {
                a = idx_dist(gen);
                b = idx_dist(gen);
                c = idx_dist(gen);
            } while (a == b || b == c || a == target);

            const size_t R = idx_dist(gen) % dim;
            for (size_t i = 0; i < dim; ++i) {
                T mut = population[a][i] + F_factor * (population[b][i] - population[c][i]);
                trials[target][i] = std::clamp(mut, lb, ub);
                if (!(uniform(gen) < CR_prob || i == R)) {
                    trials[target][i] = population[target][i];
                }
            }
        }

        for (size_t i = 0; i < population_size; ++i) {
            if (evaluation_counter >= max_evaluations) break;
            trial_fitness[i] = f(trials[i]);
            evaluation_counter++;
        }

        for (size_t i = 0; i < population_size; ++i) {
            if (trial_fitness[i] < fitness[i]) {
                population[i].swap(trials[i]);
                fitness[i] = trial_fitness[i];
                if (fitness[i] < best_fit) {
                    best_sol = population[i];
                    best_fit = fitness[i];
                }
            }
        }
    }

    return {best_sol, best_fit, evaluation_counter};
}