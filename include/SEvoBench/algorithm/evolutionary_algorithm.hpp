
#pragma once

namespace sevobench {
class evolutionary_algorithm {
protected:
  int ite = 0;
  int fes = 0;
  int mf = 0;
  int m_ite;
  int ps = 0;
  int d = 0;

public:
  evolutionary_algorithm() = default;
  evolutionary_algorithm(int _max_fes) : mf(_max_fes) {}
  evolutionary_algorithm(int _max_fes, int _pop_size)
      : mf(_max_fes), ps(_pop_size) {}
  evolutionary_algorithm(int _max_fes, int _pop_size, int _dim)
      : mf(_max_fes), ps(_pop_size), d(_dim) {}
  auto current_fes() const noexcept { return fes; }
  auto max_fes() const noexcept { return mf; }
  auto max_iterator() const noexcept { return m_ite; }
  auto set_max_iterator(int _max_iterator) noexcept { m_ite = _max_iterator; }
  auto set_max_fes(int _max_fes) noexcept { mf = _max_fes; }
  auto current_iterator() const noexcept { return ite; }
  auto dim() const noexcept { return d; }
  auto set_dim(int _dim) noexcept { d = _dim; }
  auto pop_size() const noexcept { return ps; }
  auto set_pop_size(int _pop_size) noexcept { ps = _pop_size; }
  void reset() noexcept {
    ite = 0;
    fes = 0;
  }
  auto increment_iterator() noexcept { return ++ite; }
  auto increment_fes() noexcept { return ++fes; }
  auto add_iterator(int _ite) noexcept { return ite += _ite; }
  auto add_fes(int _fes) noexcept { return fes += _fes; }
};

struct evolutionary_algorithm_condition {
  auto operator()(const auto &, const auto &alg) const noexcept {
    return alg.current_fes() <= alg.max_fes();
  }
};
} // namespace sevobench
