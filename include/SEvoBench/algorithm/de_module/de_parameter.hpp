
#pragma once

#include "../../common/tool.hpp"
#include "../evolutionary_algorithm.hpp"

namespace sevobench::de_module {
template <std::floating_point T> class de_parameter {
public:
  virtual void prepare(population<T> &) {}

  virtual void update(const population<T> &, const population<T> &) {}

  virtual T get_f(int) = 0;

  virtual T get_cr(int) = 0;

  virtual ~de_parameter() {};
};

namespace detail {
template <std::floating_point T>
auto jade_update_cr_f(const population<T> &pop, const population<T> &trial,
                      auto &&Fs, auto &&CRs, auto &&f_archive,
                      auto &&cr_archive) {
  auto f_first = f_archive.begin();
  auto cr_first = cr_archive.begin();
  auto f_begin = f_first;
  auto cr_begin = cr_first;
  for (int i = 0; i < pop.pop_size(); i++) {
    if (trial[i].fitness() < pop[i].fitness()) {
      *f_first++ = Fs[i];
      *cr_first++ = CRs[i];
    }
  }
  auto temp_c =
      (cr_begin == cr_first)
          ? T(0)
          : std::accumulate(cr_begin, cr_first, T(0)) / (cr_first - cr_begin);
  auto temp_f = (f_first == f_begin)
                    ? T(0)
                    : std::inner_product(f_begin, f_first, f_begin, T(0)) /
                          std::accumulate(f_begin, f_first, T(0));
  return std::make_pair(temp_c, temp_f);
}

template <std::floating_point T>
inline auto shade_update_cr_f(const population<T> &pop,
                              const population<T> &trial, auto &&Fs, auto &&CRs,
                              auto &&f_archive, auto &&cr_archive,
                              auto &&delta_f)

    noexcept {
  auto size = pop.pop_size();
  auto f_first = f_archive.begin();
  auto cr_first = cr_archive.begin();
  auto df_first = delta_f.begin();
  for (int i = 0; i < size; i++) {
    if (trial[i].

        fitness()

        < pop[i].

          fitness()

    ) {
      *f_first++ = Fs[i];
      *cr_first++ = CRs[i];
      *df_first++ = pop[i].

                    fitness()

                    - trial[i].

                      fitness();
    }
  }
  auto scr2 = std::inner_product(cr_archive.begin(), cr_first, delta_f.begin(),
                                 T(0), std::plus<T>(),
                                 [](auto x, auto y) { return x * x * y; });
  auto scr1 =
      std::inner_product(cr_archive.begin(), cr_first, delta_f.begin(), T(0));
  auto tmp_cr = (scr1 == T(0) ? T(-1) : scr2 / scr1);
  auto tmp_f =
      std::inner_product(f_archive.begin(), f_first, delta_f.begin(), T(0),
                         std::plus<T>(),
                         [](auto x, auto y) { return x * x * y; }) /
      std::inner_product(f_archive.begin(), f_first, delta_f.begin(), T(0));
  return std::make_pair(tmp_cr, tmp_f);
}

template <std::floating_point T>
inline void jade_prepare(const population<T> &pop, auto &&Fs, auto &&CRs,
                         auto &&RNG, T f_, T cr_)

    noexcept {
  const int size = pop.pop_size();
  auto gen_f = [&](auto x) {
    if (!std::isfinite(x))
      return T(1);
    auto z = RNG.cauchy(x, T(0.1));
    while (z <= 0)
      z = RNG.cauchy(x, T(0.1));
    return z < 1 ? z : 1;
  };
  for (int i = 0; i < size; i++) {
    CRs[i] = std::clamp(RNG.normal(cr_, T(0.1)), T(0), T(1));
    Fs[i] = gen_f(f_);
  }
}

template <std::floating_point T>
inline void shade_prepare(const population<T> &pop, auto &&Fs, auto &&CRs,
                          auto &&RNG, auto &&MF, auto &&MCR, int memory_size)

    noexcept {
  const int size = pop.pop_size();
  auto gen_f = [&](auto x) {
    if (!std::isfinite(x))
      return T(1);
    auto z = RNG.cauchy(x, T(0.1));
    constexpr int sample_max = 20;
    for (int i = 0; z <= T(0); i++) {
      z = RNG.cauchy(x, T(0.1));
      if (i > sample_max)
        return T(1);
    }
    return z < 1 ? z : T(1);
  };
  for (int i = 0; i < size; i++) {
    auto rp = RNG.rand_int(memory_size);
    CRs[i] = (MCR[rp] == T(-1)
                  ? T(0)
                  : std::clamp(RNG.normal(MCR[rp], T(0.1)), T(0), T(1)));
    Fs[i] = gen_f(MF[rp]);
  }
}

template <std::floating_point T>
inline void jso_prepare(const population<T> &pop, auto &&Fs, auto &&CRs,
                        auto &&RNG, auto &&MF, auto &&MCR, int memory_size,
                        const evolutionary_algorithm &alg)

    noexcept {
  const auto size = pop.pop_size();
  auto gen_f = [&](auto x) {
    if (!std::isfinite(x))
      return T(1);
    auto z = RNG.cauchy(x, T(0.1));
    constexpr int sample_max = 20;
    for (int i = 0; z <= T(0); i++) {
      z = RNG.cauchy(x, T(0.1));
      if (i > sample_max)
        return T(1);
    }
    return z < 1 ? z : T(1);
  };
  for (int i = 0; i < size; i++) {
    auto rp = RNG.rand_int(memory_size);
    if (rp == memory_size - 1) {
      MCR[memory_size - 1] = T(0.9);
      MF[memory_size - 1] = T(0.9);
    }
    auto tmp_cr = (MCR[rp] == T(-1)
                       ? T(0)
                       : std::clamp(RNG.normal(MCR[rp], T(0.1)), T(0), T(1)));
    if (alg.

        current_fes()

        < alg.

              max_fes()

              / 4) {
      tmp_cr = std::max(T(0.7), tmp_cr);
    } else if (alg.

               current_fes()

               < alg.

                     max_fes()

                     / 2) {
      tmp_cr = std::max(T(0.6), tmp_cr);
    }
    CRs[i] = tmp_cr;
    auto tmp_f = gen_f(MF[rp]);
    if (tmp_f > T(0.7) && alg.

                              current_fes()

                              < 3 *
                                    alg.

                                    max_fes()

                                    / 5) {
      tmp_f = T(0.7);
    }
    Fs[i] = tmp_f;
  }
}

} // namespace detail

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class jade_parameter : public de_parameter<T> {
protected:
  std::vector<T> f_archive;
  std::vector<T> cr_archive;
  std::vector<T> Fs;
  std::vector<T> CRs;
  const T c_;
  const T init_cr;
  const T init_f;
  T cr_;
  T f_;

  void resize(int size)

      noexcept {
    if (static_cast<int>(f_archive.size()) < size) {
      f_archive.resize(size);
      cr_archive.resize(size);
      Fs.resize(size);
      CRs.resize(size);
    }
  }

public:
  R RNG;

  jade_parameter(T _cr = T(0.5), T _f = T(0.5), T _c = T(0.1))
      : c_(_c), init_cr(_cr), init_f(_f), cr_(_cr), f_(_f) {}

  void prepare(population<T> &pop) override {
    resize(pop.pop_size());
    detail::jade_prepare(pop, Fs, CRs, RNG, f_, cr_);
  }

  void update(const population<T> &pop, const population<T> &trial) override {
    auto [temp_c, temp_f] =
        detail::jade_update_cr_f(pop, trial, Fs, CRs, f_archive, cr_archive);
    cr_ = (1 - c_) * cr_ + c_ * temp_c;
    f_ = (1 - c_) * f_ + c_ * temp_f;
  }

  T get_f(int i) override { return Fs[i]; }

  T get_cr(int i) override { return CRs[i]; }

  void reset()

      noexcept {
    cr_ = init_cr;
    f_ = init_f;
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class shade_parameter : public de_parameter<T> {
protected:
  const T memory_f;
  const T memory_cr;
  const int memory_size;
  std::vector<T> MF;
  std::vector<T> MCR;
  std::vector<T> f_archive;
  std::vector<T> cr_archive;
  std::vector<T> delta_f;
  std::vector<T> Fs;
  std::vector<T> CRs;
  int index_counter = 0;

  void resize(int size)

      noexcept {
    if (static_cast<int>(f_archive.size()) < size) {
      f_archive.resize(size);
      cr_archive.resize(size);
      delta_f.resize(size);
      Fs.resize(size);
      CRs.resize(size);
    }
  }

public:
  R RNG;

  shade_parameter(int _memory_size = 6, T _memory_f = T(0.5),
                  T _memory_cr = T(0.5))
      : memory_f(_memory_f), memory_cr(_memory_cr), memory_size(_memory_size),
        MF(_memory_size, _memory_f), MCR(_memory_size, _memory_cr) {}

  void prepare(population<T> &pop) override {
    resize(pop.pop_size());
    detail::shade_prepare(pop, Fs, CRs, RNG, MF, MCR, memory_size);
  }

  void update(const population<T> &pop, const population<T> &trial) override {
    auto [mcr, mf] = detail::shade_update_cr_f(pop, trial, Fs, CRs, f_archive,
                                               cr_archive, delta_f);
    MCR[index_counter] = mcr;
    MF[index_counter] = mf;
    index_counter = (index_counter + 1) % memory_size;
  }

  T get_f(int i) override { return Fs[i]; }

  T get_cr(int i) override { return CRs[i]; }

  void reset()

      noexcept {
    index_counter = 0;
    std::fill(MCR.begin(), MCR.end(), memory_cr);
    std::fill(MF.begin(), MF.end(), memory_f);
  }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class jso_parameter final : public shade_parameter<T, R> {
  const evolutionary_algorithm &alg;

public:
  jso_parameter(const evolutionary_algorithm &_alg, int _memory_size = 5,
                T _memory_f = T(0.3), T _memory_cr = T(0.8))
      : shade_parameter<T, R>(_memory_size, _memory_f, _memory_cr), alg(_alg) {}

  void prepare(population<T> &pop) override {
    this->resize(pop.pop_size());
    detail::jso_prepare(pop, this->Fs, this->CRs, this->RNG, this->MF,
                        this->MCR, this->memory_size, alg);
  }
};

template <std::floating_point T>
class constant_parameter final : public de_parameter<T> {
  const T f;
  const T cr;

public:
  constant_parameter(T _f = T(0.5), T _cr = T(0.9)) : f(_f), cr(_cr) {}

  T get_f(int) override { return f; }

  T get_cr(int) override { return cr; }
};

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>

class jde_parameter : public de_parameter<T> {
  const T init_f;
  const T init_cr;
  const T fl;
  const T fu;
  const T t1;
  const T t2;
  T f;
  T cr;

public:
  R RNG1;
  R RNG2;
  R RNG3;
  R RNG4;

  jde_parameter(T _fl = T(0.1), T _fu = T(0.9), T _t1 = T(0.1), T _t2 = T(0.1),
                T _init_f = T(0.5), T _init_cr = T(0.9))
      : init_f(_init_f), init_cr(_init_cr), fl(_fl), fu(_fu), t1(_t1), t2(_t2) {
  }

  T get_f(int) override {
    f = RNG2.template rand_float<T>() < t1
            ? fl + fu * RNG1.template rand_float<T>()
            : f;
    return f;
  }

  T get_cr(int) override {
    cr =
        RNG4.template rand_float<T>() < t2 ? RNG3.template rand_float<T>() : cr;
    return cr;
  }

  void reset()

      noexcept {
    f = init_f;
    cr = init_cr;
  }
};

} // namespace sevobench::de_module
