#pragma once
#include "eoRNG.h"
class paradiseo_rng {
    mutable eoRng sr;

public:
    paradiseo_rng(std::uint32_t _ =
                          [] {
                              std::random_device rd;
                              return rd();
                          }())
        : sr(_) {};
    template<typename T>
    T rand_float(T l = T(0), T r = T(1)) const noexcept {
        return sr.uniform(l, r);
    }

    auto operator()() const noexcept { return sr.rand(); }

    int rand_int(int n) const noexcept { return sr.random(n); }
    template<typename T>
    T normal(T m, T st) const noexcept {
        return box_muller(m, st, rand_float<T>(), rand_float<T>());
    }
    template<typename T>
    T cauchy(T a, T b) const noexcept {
        return cauchy_dis(a, b, rand_float<T>());
    }
    void seed(unsigned int _) const noexcept { sr.reseed(_); }

    template<int K>
    auto pick_random(int n, int j) const noexcept {
        std::array<int, K> select{};
        for (int i = 0; i < K; i++) {
            do {
                select[i] = rand_int(n);
            } while (select[i] == j ||
                     std::any_of(select.begin(), select.begin() + i,
                                 [&](auto x) { return select[i] == x; }));
        }
        return select;
    }
};