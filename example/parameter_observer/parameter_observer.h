#pragma once
#include "SEvoBench/sevobench.hpp"
namespace sevobench::de_module {

template <std::floating_point T, typename R = tool::rng>
  requires tool::random_generator_concept<R, T>
class trackable_jade_parameter : public jade_parameter<T, R> {
public:
  using jade_parameter<T, R>::jade_parameter; // Inherit constructors

  //  get_f() to return the internal f_ value
  T get_f() { return this->f_; }

  // get_cr() to return the internal cr_ value
  T get_cr() { return this->cr_; }
};

} // namespace sevobench::de_module

namespace sevobench::experiment {

template <std::floating_point T>
class jade_parameter_observer final : public suite_observer<T> {
public:
  using value_type = T;
  using parameter_data_type =
      std::vector<std::array<T, 2>>; // [f_value, cr_value]

private:
  std::vector<std::vector<parameter_data_type>>
      data_; // Stores parameter changes for each problem, instance, and run
  const std::vector<int> pro_index; // Problem indices in the suite
  const int ins_count_;             // Number of instances per problem
  const int max_fes_;               // Maximum number of fitness evaluations
  const int runs_;                  // Number of independent runs
  de_module::trackable_jade_parameter<T>
      *jade_param_;     // Pointer to the trackable JADE parameter object
  int population_size_; // Population size (fixed during the optimization
                        // process)

public:
  jade_parameter_observer(const auto &_suite, int _max_fes, int _runs,
                          de_module::trackable_jade_parameter<T> *jade_param,
                          int population_size)
      : pro_index(_suite.problem_index()), ins_count_(_suite.instance_count()),
        max_fes_(_max_fes), runs_(_runs), jade_param_(jade_param),
        population_size_(population_size) {
    // Resize data_ to match the number of problems and instances
    data_.resize(_suite.size());
    for (auto &d : data_) {
      d.resize(_runs);
      for (auto &v : d) {
        // Reserve space based on the maximum number of iterations
        v.reserve(_max_fes / population_size_);
      }
    }
  }

  void log(const problem::problem_state<T> &state,
           const problem::problem_info<T> &info) override {
    // Calculate the problem index (m) and instance index (n)
    auto m = std::find(pro_index.begin(), pro_index.end(), info.index) -
             pro_index.begin();
    auto n = info.instance;
    // auto p = (m * ins_count_ + n - 1) * runs_ + state.run_id - 1;

    // Check if this is the start of a new iteration
    if (state.evaluations % population_size_ == 0 && state.evaluations > 0) {
      if (state.evaluations <= max_fes_) {
        auto i1 = (m * ins_count_ + n - 1);
        T f_value = jade_param_->get_f();   // Get the current f_ value
        T cr_value = jade_param_->get_cr(); // Get the current cr_ value
        // Record f_value and cr_value
        data_[i1][state.run_id - 1].push_back({f_value, cr_value});
      }
    }
  }

  // Get parameter data for a specific problem index
  auto at(int problem_index) const noexcept {
    auto m = std::find(pro_index.begin(), pro_index.end(), problem_index) -
             pro_index.begin();
    return std::span<const std::vector<parameter_data_type>>(
        data_.data() + m * ins_count_, ins_count_);
  }

  // Get parameter data for a specific problem index and instance
  auto at(int problem_index, int instance) const noexcept {
    auto m = std::find(pro_index.begin(), pro_index.end(), problem_index) -
             pro_index.begin();
    return std::span<const parameter_data_type>(
        data_[m * ins_count_ + instance - 1]);
  }
};

} // namespace sevobench::experiment
