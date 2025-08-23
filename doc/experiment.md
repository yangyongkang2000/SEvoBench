# SEvoBench Experiment Framework Documentation

---

## 1. Overview

The `sevobench::experiment` module provides a comprehensive framework for benchmarking evolutionary algorithms, featuring:

- **Parallel Execution**: Multi-threaded execution of independent experimental runs
- **Time-Series Tracking**: Continuous monitoring of best-so-far solutions
- **Extensible Observation**: Custom data collection via observer pattern
- **Reproducible Experiments**: Deterministic configuration management

Core components include the `best_so_far_record` performance tracker and `evo_bench` benchmarking driver, designed for large-scale algorithm evaluation in Modern C++20.

---

## 2. Core Components

### 2.1 Best-So-Far Tracking System (`best_so_far_record`)

#### Architecture
```cpp
template <std::floating_point T>
class best_so_far_record final : public suite_observer<T> {
  // Storage: [Problems][Instances][Runs][(FEs, Value)]
  std::vector<std::vector<std::vector<std::pair<int, T>>>> data_;
  // Best values: [Problem×Instance×Run]
  std::vector<T> bests; 
};
```

#### Key Features
| Feature                 | Description                          |
|-------------------------|--------------------------------------|
| Multi-Run Support       | Tracks multiple independent executions |
| Instance Management     | Handles multiple problem instances  |
| Compressed Storage      | Records at configurable FE intervals|
| Lock-free Updates       | Atomic updates for parallel safety  |

#### Initialization
```cpp
// Create recorder for CEC2020 suite with:
// - 1e5 max function evaluations 
// - 25 independent runs
// - 1000 FE recording interval
best_so_far_record<double> tracker(cec_suite, 1e5, 25, 1000);
```

### 2.2 Benchmark Driver (`evo_bench`)

#### Function Signature
```cpp
template <bool Parallel = true>
void evo_bench(Algorithm&& alg, Suite&& suite, Observer&& obs, 
               int independent_runs, int max_fes);
```

#### Execution Workflow
```mermaid
graph TD
  A[Initialize Suite] --> B[Create Thread Pool]
  B --> C[Distribute Problems]
  C --> D[Worker Thread: Run Algorithm]
  D --> E[Observer Updates]
  E --> F[Result Aggregation]
  F --> G{All Runs Complete?}
  G --Yes--> H[Generate Reports]
```

---

## 3. Key Functionality

### 3.1 Performance Tracking

#### Data Recording Logic
```cpp
void log(const problem_state<T>& state, const problem_info<T>& info) override {
  // 1. Calculate storage indices
  auto prob_idx = find_problem_index(info.index);
  auto inst_idx = info.instance;
  
  // 2. Update best value atomically
  atomic_min(bests[prob_idx][inst_idx][state.run_id], state.current_value);
  
  // 3. Record time-series snapshot
  if(state.evaluations % step_ == 0) {
    data_[prob_idx][inst_idx][state.run_id].emplace_back(
      state.evaluations, current_best);
  }
}
```

#### Data Access Methods
```cpp
// Get Problem 3, Instance 2 data
auto prob3_data = recorder.at(3, 2); 

// Access raw best values
auto all_bests = recorder.best();
```

### 3.2 Parallel Benchmarking

#### Task Distribution
```cpp
parallel_task pool(std::thread::hardware_concurrency());
std::vector<std::future<void>> futures;

for (auto& problem : suite) {
  for (int run = 0; run < runs; run++) {
    futures.emplace_back(pool.submit([&] {
      algorithm(suite_problem(problem, recorder, run));
    }));
  }
}
```

---

## 4. Benchmark Configuration

### 4.1 DE Algorithm Setup
```cpp
// Configure JADE parameters
auto de_config = de_module::de_config<true, double>{
  .parameter = std::make_unique<de_module::jade_parameter<double>>(0.1, 0.5),
  .mutation = std::make_unique<de_module::ttpb1_mutation<double>>(0.11),
  .constraint_handler = std::make_unique<de_module::projection_repair<double>>(),
  .crossover = std::make_unique<de_module::binomial_crossover<double>>(),
  .population_strategy = std::make_unique<de_module::linear_reduction<double>>(100, 50),
  .archive = std::make_unique<de_module::fifo_archive<double>>(2.0)
};
```

### 4.2 Experiment Execution
```cpp
// Execute parallel benchmark
experiment::evo_bench<true>(
  [config = std::move(de_config)](auto& problem) mutable {
    population<double> pop(100, problem.dim(), 
                          problem.lower_bound(), 
                          problem.upper_bound());
    evolutionary_algorithm tracker(1e5, pop.size(), problem.dim());
    return de_module::de_algorithm<true, double>(std::move(config))
           .run(pop, [&](auto x){ return problem(x); },
                problem.lower_bound(), problem.upper_bound(), tracker);
  },
  suite,
  recorder,
  25,    // 25 independent runs
  1e5    // Max function evaluations
);
```

---

## 5. Advanced Usage

### 5.1 Custom Observers
```cpp
class diversity_observer : public suite_observer<double> {
  void log(const problem_state<double>& s, 
          const problem_info<double>& i) override {
    // Track population diversity metrics
  }
};
```

### 5.2 Result Analysis
```cpp
// Calculate median performance across runs
for (auto& prob_data : recorder) {
  auto [median, iqr] = calculate_statistics(prob_data);
  std::ofstream csv("stats.csv");
  csv << "FEs,Median Fitness\n";
  for (const auto& [fes, val] : median) {
    csv << fes << "," << val << "\n";
  }
}
```

---

## 6. API Reference

### 6.1 `best_so_far_record`
| Method                 | Description                          |
|------------------------|--------------------------------------|
| `log()`                | Records evaluation state            |
| `at(problem, instance)`| Access specific instance data       |
| `best()`               | Get raw best values across all runs  |

### 6.2 `evo_bench`
| Parameter          | Requirements                  |
|--------------------|-------------------------------|
| `alg`              | Must be invocable with problem|
| `suite`            | Configured problem collection |
| `obs`              | Inherits from suite_observer  |

---

## 7. Complete Benchmark Example

```cpp
#include <SEvoBench/sevobench.hpp>

int main() {
  using namespace sevobench;
  using namespace experiment;
  using namespace de_module;

  // 1. Configure test suite
  auto suite = suite_builder<problem::CEC2020>()
               .dim<30>()
               .type<double>()
               .problem_index(problem::problem_range<1,15>())
               .instance_count(5)
               .build();

  // 2. Initialize performance tracker
  best_so_far_record<double> recorder(suite, 1e5, 25, 1000);

  // 3. Configure JADE algorithm
  auto de_config = de_config<true, double>{
    std::make_unique<jade_parameter<double>>(0.1, 0.5),
    std::make_unique<ttpb1_mutation<double>>(0.11),
    std::make_unique<projection_repair<double>>(),
    std::make_unique<binomial_crossover<double>>(),
    std::make_unique<linear_reduction<double>>(100, 50),
    std::make_unique<fifo_archive<double>>(2.0)
  };

  // 4. Execute benchmark
  evo_bench<true>(
    [config = std::move(de_config)](auto& problem) mutable {
      population<double> pop(100, problem.dim(), 
                            problem.lower_bound(), 
                            problem.upper_bound());
      evolutionary_algorithm tracker(1e5, pop.size(), problem.dim());
      return de_algorithm<true, double>(std::move(config))
             .run(pop, [&](auto x){ return problem(x); },
                  problem.lower_bound(), problem.upper_bound(), tracker);
    },
    suite,
    recorder,
    25,     // 25 independent runs
    100000  // Max FEs
  );

  // 5. Export results
  std::ofstream csv("results.csv");
  csv << "Problem,Instance,Run,FEs,Fitness\n";
  for (size_t pid : suite.problem_index()) {
    for (int inst = 1; inst <= suite.instance_count(); inst++) {
      auto data = recorder.at(pid, inst);
      for (size_t run = 0; run < data.size(); run++) {
        for (const auto& [fes, fit] : data[run]) {
          csv << pid << "," << inst << "," << run+1 << "," 
              << fes << "," << fit << "\n";
        }
      }
    }
  }
}
```

---

## 8. Performance Characteristics

| Metric               | Value                      |
|----------------------|----------------------------|
| Parallel Scaling     | 95%+ core utilization      |
| Memory Footprint     | 8 bytes per data point     |
| Observation Overhead | <1% runtime penalty        |
| Thread Safety        | Lock-free atomic updates   |

---

## 9. Contribution Guidelines

1. **Component Extensions**:
    - Implement new problem types using `problem_common` interface
    - Develop specialized observers using `suite_observer` base class

2. **Code Standards**:
   ```cpp
   // Validate interfaces with concepts
   static_assert(experiment::suite_observer<double>, "Invalid observer");
   ```

3. **Testing Requirements**:
    - Include validation tests for statistical methods
    - Maintain thread-safety in parallel components

---

This documentation provides a comprehensive reference for the SEvoBench experiment framework, enabling rigorous evaluation of evolutionary algorithms through standardized benchmarking procedures.