#pragma once
#include "experiment.hpp"
namespace sevobench::experiment {
template <std::floating_point T>
class best_so_far_record final : public suite_observer<T> {
public:
  using value_type = T;
  using data_type = std::vector<std::pair<int, T>>;

private:
  std::vector<T> bests;
  const std::vector<int> pro_index;
  std::vector<std::vector<data_type>> data_;
  const int ins_count_;
  const int max_fes_;
  const int step_;
  const int runs_;

public:
  best_so_far_record(const auto &_suite, int _max_fes, int _runs, int _step)
      : bests(_runs * _suite.size(), std::numeric_limits<T>::max()),
        pro_index(_suite.problem_index()), ins_count_(_suite.instance_count()),
        max_fes_(_max_fes), step_(_step), runs_(_runs) {
    data_.resize(_suite.size());
    for (auto &d : data_) {
      d.resize(_runs);
      for (auto &v : d) {
        v.reserve(_max_fes / _step);
      }
    }
  }
  void log(const problem::problem_state<T> &state,
           const problem::problem_info<T> &info) override {
    auto m = std::find(pro_index.begin(), pro_index.end(), info.index) -
             pro_index.begin();
    auto n = info.instance;
    auto p = (m * ins_count_ + n - 1) * runs_ + state.run_id - 1;
    auto &_best = bests[p];
    _best = std::min(_best, state.current_value);
    if (state.evaluations % step_ == 0) {
      if (state.evaluations <= max_fes_) {
        auto i1 = (m * ins_count_ + n - 1);
        data_[i1][state.run_id - 1].emplace_back(state.evaluations, _best);
      }
    }
  }
  auto at(int problem_index) const noexcept {
    auto m = std::find(pro_index.begin(), pro_index.end(), problem_index) -
             pro_index.begin();
    return std::span<const std::vector<data_type>>(
        data_.data() + m * ins_count_, ins_count_);
  }
  auto at(int problem_index, int instance) const noexcept {
    auto m = std::find(pro_index.begin(), pro_index.end(), problem_index) -
             pro_index.begin();
    return std::span<const data_type>(data_[m * ins_count_ + instance - 1]);
  }
  auto best() const noexcept { return std::span<const T>(bests); }
};
} // namespace sevobench::experiment