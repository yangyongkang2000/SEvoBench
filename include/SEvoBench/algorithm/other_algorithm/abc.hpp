#pragma once
#include "../../common/common_concept.hpp"
#include "../../common/tool.hpp"
namespace sevobench::other_algorithm::abc_detail {

    template<int Dim, int Pop_Size, std::floating_point T>
    struct colony {
        using value_type = T;
        static constexpr int m = Pop_Size;
        const int l;
        static constexpr int n = Dim;
        static constexpr T qk = T(1) / T(0x7FFF);
        int t[Pop_Size]{}, tmp[3 * Pop_Size]{};
        T fit[2 * Pop_Size]{};
        tool::simple_rand rd;
        T d, w;
        colony(unsigned int x, T _d, T _w, int _l)
            : l(_l), rd(tool::simple_rand(x)), d(_d), w(_w) {}
    };
    template<typename F, typename InputIt>
    inline auto function_fitness(F &&f, InputIt beg) noexcept {
        auto tmp = f(beg);
        return tmp >= 0 ? 1 / (1 + tmp) : 1 - tmp;
    }
    template<typename S, typename F, typename C>
    inline void colony_init(S &sol, F &&f, C &c) noexcept {
        std::transform(sol.begin(), sol.end(), c.fit,
                       [&](auto &v) { return function_fitness(f, v.data()); });
    }
    template<typename S, typename F, typename C>
    inline void colony_search(S &sol, F &&f, C &c, int m, int k, int i,
                              int q) noexcept {
        auto xmi = sol[m][i];
        auto xki = sol[k][i];
        sol[m][i] += (2 * q * C::qk - 1) * (xmi - xki);
        if (sol[m][i] > c.w)
            sol[m][i] = c.w;
        if (sol[m][i] < c.d)
            sol[m][i] = c.d;
        auto v = function_fitness(f, sol[m].data());
        if (v < c.fit[m]) {
            c.t[m]++;
            sol[m][i] = xmi;
        } else {
            c.fit[m] = v;
            c.t[m] = 0;
        }
    }
    template<typename S, typename F, typename C>
    inline void colony_employ(S &sol, F &&f, C &c) noexcept {
        std::generate_n(c.tmp, 3 * C::m, [&] { return c.rd(); });
        for (int m = 0; m < C::m; m++)
            colony_search(sol, f, c, m,
                          (c.tmp[3 * m] % C::m) == m ? (m + 1) % C::m
                                                     : c.tmp[3 * m] % C::m,
                          c.tmp[3 * m + 1] % C::n, c.tmp[3 * m + 2]);
    }
    template<typename S, typename F, typename C>
    inline void colony_outlook(S &sol, F &&f, C &c) noexcept {
        std::inclusive_scan(c.fit, c.fit + C::m, c.fit + C::m);
        auto t = c.rd() * C::qk * c.fit[2 * C::m - 1];
        auto m = static_cast<int>(
                std::lower_bound(c.fit + C::m, c.fit + 2 * C::m, t) - c.fit - C::m);
        auto k = c.rd() % C::m;
        colony_search(sol, f, c, m, k == m ? (k + 1) % C::m : k, c.rd() % C::n,
                      c.rd());
    }
    template<typename S, typename F, typename C>
    inline int colony_scout(S &sol, F &&f, C &c) noexcept {
        auto it = std::max_element(c.t, c.t + C::m);
        if (*it > c.l) {
            *it = 0;
            auto index = it - c.t;
            std::generate_n(sol[index].begin(), C::n,
                            [&] { return c.d + (c.w - c.d) * C::qk * c.rd(); });
            c.fit[index] = function_fitness(f, sol[index].data());
            return 1;
        }
        return 0;
    }
    template<typename S, typename C, typename R>
    inline void max_colony(S &sol, C &c, R &r) noexcept {
        auto it = std::max_element(c.fit, c.fit + C::m);
        auto index = it - c.fit;
        if (r.second < *it) {
            r.second = *it;
            std::copy_n(sol[index].begin(), C::n, r.first.begin());
        }
    }

}// namespace sevobench::other_algorithm::abc_detail

namespace sevobench::other_algorithm {

    template<std::floating_point T>
    struct abc_parameter {
        T ratio = T(0.5);
        std::optional<int> max{};
    };

    template<typename P, typename T>
    concept abc_parameter_concept = requires(const P &p) {
        requires std::same_as<std::remove_cvref_t<decltype(p.ratio)>, T>;
        requires std::same_as<std::optional<int>, std::remove_cvref_t<decltype(p.max)>>;
    };

    template<int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
             bool Memory_Flag = false, typename G, typename F,
             std::floating_point T, typename Parameter_Type = abc_parameter<T>>
        requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
                 abc_parameter_concept<Parameter_Type, T> &&
                 algorithm_positions_concept<F, G, T>
    inline auto abc_optimize(G &&sol, F &&f, T d, T w,
                             const Parameter_Type &pt = Parameter_Type()) noexcept {
        tool::curve_t<Memory_Flag, T, Max> abc_convergence_curve;
        std::random_device rd;
        std::pair<std::array<T, Dim>, T> result{{}, {}};
        abc_detail::colony<Dim, Pop_Size, T> c(
                rd(), d, w, static_cast<int>(pt.ratio * Dim * Pop_Size));
        abc_detail::colony_init(sol, f, c);
        int fes(Pop_Size);
        const int REAL_MAX = pt.max ? pt.max.value() : Max;
        for (int i = 0; (Memory_Flag ? i : fes) < REAL_MAX; i++) {
            abc_detail::colony_employ(sol, f, c);
            abc_detail::colony_outlook(sol, f, c);
            abc_detail::max_colony(sol, c, result);
            fes += (abc_detail::colony_scout(sol, f, c) + Pop_Size + 1);
            if constexpr (Memory_Flag) {
                abc_convergence_curve[2 * i] = T(fes);
                abc_convergence_curve[2 * i + 1] =
                        result.second <= T(1) ? T(1) / result.second - 1 : 1 - result.second;
            }
        }
        if constexpr (Memory_Flag)
            return std::make_tuple(result.first,
                                   result.second <= T(1) ? T(1) / result.second - 1
                                                         : 1 - result.second,
                                   std::move(abc_convergence_curve));
        if constexpr (!Memory_Flag)
            return std::make_pair(result.first, result.second <= T(1)
                                                        ? T(1) / result.second - 1
                                                        : 1 - result.second);
    }

    template<int Dim, int Pop_Size = 100, int Max = 1000 * Dim,
             bool Memory_Flag = false, typename F, std::floating_point T,
             typename Parameter_Type = abc_parameter<T>>
        requires algorithm_func_concept<Dim, Pop_Size, Max, F, T> &&
                 abc_parameter_concept<Parameter_Type, T>
    inline auto abc(F &&f, T d, T w,
                    const Parameter_Type &pt = Parameter_Type()) noexcept {
        auto pop(tool::random_generate_position<T, Dim, Pop_Size>(d, w));
        return abc_optimize<Dim, Pop_Size, Max, Memory_Flag>(pop, f, d, w, pt);
    }

}// namespace sevobench::other_algorithm
