//
// Created by 杨永康 on 2024/2/14.
//

#ifndef ADD_PROBLEM_ADD_PROBLEM_HPP
#define ADD_PROBLEM_ADD_PROBLEM_HPP

#include "SEvoBench/single_problem.hpp"

template <int Prob_Index, int Dim, typename T> class add_problem;

template <int Dim, typename T>
  requires(Dim >= 4 && Dim % 2 == 0)
class add_problem<1, Dim, T> {
public:
  static constexpr T L = 0;
  static constexpr T U = 1;

  T operator()(T *points) {
    constexpr int N = Dim / 2;
    T minDistance = std::numeric_limits<T>::infinity();
    for (int i = 0; i < N - 1; ++i) {
      T xi = points[2 * i];
      T yi = points[2 * i + 1];
      for (int j = i + 1; j < N; ++j) {
        T xj = points[2 * j];
        T yj = points[2 * j + 1];
        minDistance = std::min((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi),
                               minDistance);
      }
    }
    return 2 - minDistance;
  }
};

template <typename T, int N>
void init(){

};

REGISTER_PROBLEMS(AddProblem, add_problem, 1)

INIT_PROBLEMS(add_problem, init)




#endif // ADD_PROBLEM_ADD_PROBLEM_HPP
