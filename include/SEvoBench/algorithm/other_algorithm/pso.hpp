#pragma once
#include "../../common/common_concept.hpp"
#include "../../common/tool.hpp"

namespace sevobench::other_algorithm {

    template<std::floating_point T>
    struct pso_parameter {
        T c1 = T(2);
        T c2 = T(2);
        T w_min = T(0.4);
        T w_max = T(0.9);
        T v_max_ratio = T(0.2);
        std::optional<int> max{};
    };

    template<typename P, typename T>
    concept pso_parameter_concept = requires(const P &p) {
        requires std::same_as<T, std::remove_cvref_t<decltype(p.c1)>>;
        requires std::same_as<T, std::remove_cvref_t<decltype(p.c2)>>;
        requires std::same_as<T, std::remove_cvref_t<decltype(p.w_min)>>;
        requires std::same_as<T, std::remove_cvref_t<decltype(p.w_max)>>;
        requires std::same_as<T, std::remove_cvref_t<decltype(p.v_max_ratio)>>;
        requires std::same_as<std::optional<int>, std::remove_cvref_t<decltype(p.max)>>;
    };

    template<int Dim, int Pop_Size = 30, int Max = 1000 * Dim,
             bool Memory_Flag = false, typename G, typename F,
             std::floating_point T, typename Parameter_Type = pso_parameter<T>>
        requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
                 pso_parameter_concept<Parameter_Type, T> &&
                 algorithm_positions_concept<F, G, T>
    inline auto pso_optimize(G &&positions, F &&f, T left_bound, T right_bound,
                             const Parameter_Type &pt = Parameter_Type()) noexcept {
        const T k = T(1) / T(0x7FFF);
        const T c1 = pt.c1;
        const T c2 = pt.c2;
        const T w_max = pt.w_max;
        const T w_min = pt.w_min;
        T v_max = (right_bound - left_bound) * pt.v_max_ratio;
        std::random_device rd;
        std::default_random_engine e(rd());
        tool::simple_rand sr1(rd());
        tool::simple_rand sr2(rd());
        std::uniform_real_distribution<T> v(-v_max, v_max);
        std::array<T, Dim> r1;
        std::array<T, Dim> r2;
        tool::curve_t<Memory_Flag, T, Max> pso_convergence_curve;
        auto particles_best_pos(positions);
        std::array<T, Pop_Size> particles_best_score;
        std::transform(positions.begin(), positions.end(),
                       particles_best_score.begin(),
                       [&](auto &x) { return f(x.data()); });
        auto global_best_index = std::min_element(particles_best_score.begin(),
                                                  particles_best_score.end()) -
                                 particles_best_score.begin();
        auto global_best_score = particles_best_score[global_best_index];
        std::vector<std::array<T, Dim>> velocities(Pop_Size);
        for (auto &x: velocities)
            std::generate_n(x.begin(), Dim, [&] { return v(e); });
        const int REAL_MAX = pt.max ? pt.max.value() : Max;
        for (int iteration{Memory_Flag ? 0 : Pop_Size};
             iteration < REAL_MAX;
             iteration += (Memory_Flag ? 1 : Pop_Size)) {
            auto w = w_max - iteration * ((w_max - w_min) / REAL_MAX);
            for (int i = 0; i < Pop_Size; i++) {
                for (int j = 0; j < Dim; j++) {
                    r1[j] = sr1() * k;
                    r2[j] = sr2() * k;
                }
                for (int j = 0; j < Dim; j++)
                    velocities[i][j] = std::clamp(
                            w * velocities[i][j] +
                                    c1 * r1[j] * (particles_best_pos[i][j] - positions[i][j]) +
                                    c2 * r2[j] *
                                            (particles_best_pos[global_best_index][j] -
                                             positions[i][j]),
                            -v_max, v_max);
            }
            for (int i = 0; i < Pop_Size; i++) {
                for (int j = 0; j < Dim; j++)
                    positions[i][j] = std::clamp(positions[i][j] + velocities[i][j],
                                                 left_bound, right_bound);
                auto tmp = f(positions[i].data());
                if (tmp < particles_best_score[i]) {
                    particles_best_pos[i] = positions[i];
                    particles_best_score[i] = tmp;
                }
            }
            global_best_index = std::min_element(particles_best_score.begin(),
                                                 particles_best_score.end()) -
                                particles_best_score.begin();
            global_best_score = particles_best_score[global_best_index];
            if constexpr (Memory_Flag) {
                pso_convergence_curve[2 * iteration] = (iteration + 2) * Pop_Size;
                pso_convergence_curve[2 * iteration + 1] = global_best_score;
            }
        }
        if constexpr (Memory_Flag)
            return std::make_tuple(positions[global_best_index], global_best_score,
                                   std::move(pso_convergence_curve));
        if constexpr (!Memory_Flag)
            return std::make_pair(positions[global_best_index], global_best_score);
    }

    template<int Dim, int Pop_Size = 30, int Max = 1000 * Dim,
             bool Memory_Flag = false, typename F, std::floating_point T,
             typename Parameter_Type = pso_parameter<T>>
        requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
                 pso_parameter_concept<Parameter_Type, T>
    inline auto pso(F &&f, T l, T r,
                    const Parameter_Type &pt = Parameter_Type()) noexcept {
        auto pop(tool::random_generate_position<T, Dim, Pop_Size>(l, r));
        return pso_optimize<Dim, Pop_Size, Max, Memory_Flag>(pop, f, l, r, pt);
    }

}// namespace sevobench::other_algorithm
