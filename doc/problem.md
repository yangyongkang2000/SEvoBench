# SEvoBench Problem Module Documentation

The `sevobench::problem` module provides a unified interface for defining, configuring, and managing optimization problems within the SEvoBench framework. This document outlines its architecture, problem management patterns, and core capabilities.

---

## 1. Overview
The Problem module implements:
- **Problem Suite Management**: Unified container for optimization problems
- **Generic Problem Interface**: Standardized evaluation interface
- **Configuration System**: Type-safe problem suite construction
- **Instance Generation**: Official/randomized problem instances
- **Meta Information**: Bounds, optima, and problem statistics
- **Dimensionality Support**: Flexible dimension configuration

Key features include multi-problem benchmarking, instance randomization, and automatic data loading for reproducible research.

---

## 2. Core Components

### 2.1 Problem Base Classes
| Class                    | Purpose                          |
|--------------------------|----------------------------------|
| `problem_common`         | CRTP base for problem definitions |
| `single_problem`         | Abstract interface for problems  |
| `problem_info`           | Metadata storage (bounds/optima) |

### 2.2 Suite Management
| Component               | Functionality                   |
|-------------------------|---------------------------------|
| `suite`                 | Problem collection container    |
| `suite_builder`         | Type-safe configuration builder |
| `problem_range`         | Generate problem ID sequences   |

---

## 3. Configuration System

### 3.1 Suite Builder Pattern
```cpp
auto builder = problem::suite_builder<ProblemType>()
                 .dim<30>()               // 30-dimensional problems
                 .type<double>()         // Double-precision evaluation
                 .problem_index({1,3,5}) // Select problems F1/F3/F5
                 .instance_count(5)       // Generate 5 randomized instances
                 .dir("official_data");   // Load official instance data

auto problem_suite = builder.build();
```

### 3.2 Configuration Constraints
| Method                  | Requirement                  | Error Condition              |
|-------------------------|------------------------------|------------------------------|
| `.dim<N>()`             | Must specify dimensions      | Compile-time assertion       |
| `.type<T>()`            | Must specify floating type   | Compile-time type check      |
| `.problem_index()`      | Valid problem IDs required   | Runtime clamping/sorting     |
| `.instance_count()`     | ≥1 instances                 | Auto-correct to 1 if invalid |

---

## 4. Key Features

### 4.1 Problem Interface
**Standard Evaluation Protocol:**
```cpp
T operator()(std::span<const T> x) {
  // Problem-specific evaluation logic
  return fitness;
}
```

**Metadata Access:**
```cpp
auto info = problem.problem_information();
std::cout << "Bounds: [" << info.lb << ", " << info.ub << "]\n"
          << "Optimum: " << info.optimum.value() << "\n";
```

### 4.2 Instance Management
**Official Instances:**
```cpp
suite<MyProblem> official_suite({1,2,3}, "data/"); // Load preconfigured data
```

**Random Instances:**
```cpp
suite<MyProblem> random_suite({4,5}, 10); // 10 randomized instances
```

### 4.3 Dimensionality Handling
**Fixed-Dimension Support:**
```cpp
template <int Dim, typename T>
class MyProblem : public problem_common<1, Dim, T, MyProblem> { ... };
```

**Dimension Validation:**
```cpp
static_assert(Dim % 10 == 0, "Dimension must be multiple of 10");
```

---

## 5. API Reference

### 5.1 Core Classes
| Class                    | Key Methods                     |
|--------------------------|---------------------------------|
| `suite`                  | `begin()`, `end()`, `size()`    |
| `suite_builder`          | `dim()`, `type()`, `build()`     |
| `single_problem`         | `operator()`, `optimum()`        |

### 5.2 Key Interfaces
| Method                          | Description                          |
|---------------------------------|--------------------------------------|
| `problem_information()`         | Returns bounds/optima metadata       |
| `problem_range<First,Last>()`   | Generates ID sequence [First,Last]  |
| `generate_problem_factory()`    | Creates problem instance generators  |

---

## 6. Usage Example

```cpp
#include "SEvoBench/sevobench.hpp"

int main() {
  using namespace sevobench::problem;
  
  // Configure benchmark suite
  auto suite = suite_builder<MyProblem>()
                 .dim<50>()
                 .type<float>()
                 .problem_index({2,4,6})
                 .instance_count(3)
                 .build();

  // Evaluate solutions
  std::vector<float> x(50, 0.5f);
  for (auto& prob : suite) {
    auto fitness = prob(x);
    std::cout << "F" << prob.index() 
              << " Instance " << prob.instance()
              << " Fitness: " << fitness << "\n";
  }
}
```

---

## 7. Advanced Features

### 7.1 Problem Factories
**Generator Creation:**
```cpp
auto factory = generate_problem_factory<MyProblem, 30, double>();
auto prob = factory[0](1); // Create first problem instance
```

**Custom Initialization:**
```cpp
template <int Dim, typename T>
class CustomProblem : public problem_common<7, Dim, T, CustomProblem> {
  CustomProblem(int instance) {
    // Custom instance initialization
  }
};
```

### 7.2 Dynamic Configuration
**Runtime Problem Selection:**
```cpp
std::vector<int> active_problems;
// ... load from config file
auto suite = suite_builder<MyProblem>()
               .problem_index(active_problems)
               .build();
```

---

## 8. Performance Considerations

- **Evaluation Overhead**: Virtual dispatch penalty <1% for problem dimensions >10
- **Memory Footprint**:
    - Base problem: 24 bytes + metadata
    - 1000-problem suite: ~50KB overhead
- **Parallel Evaluation**: Thread-safe for concurrent fitness evaluation

---

## 9. Contribution Guidelines

1. Implement new problems via `problem_common` CRTP pattern
2. Maintain standardized evaluation interface:
   ```cpp
   T operator()(std::span<const T> x) const noexcept;
   ```
3. Provide metadata through `problem_info` structure
4. Support both fixed and configurable dimensions
5. Include instance generation logic:
    - Official data loading
    - Randomized initialization
6. Validate against reference implementations
7. Document problem characteristics in header comments

---

This modular design enables researchers to integrate new optimization problems while maintaining compatibility with SEvoBench's benchmarking infrastructure. The type-safe builder pattern ensures correct configurations while preserving flexibility for experimental setups.