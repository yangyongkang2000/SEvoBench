#ifndef problem_wrapper_h
#define problem_wrapper_h

extern "C" {
#include "coco.h"
}
#include "SEvoBench/single_algorithm.hpp"
#include <cstdio>

class problem_wrapper {
  coco_problem_t *problem;

public:
  problem_wrapper(coco_problem_t *_p) : problem(_p){};
  double operator()(double *x) {
    double y;
    coco_evaluate_function(problem, x, &y);
    return y;
  }
  double lower() {
    return *coco_problem_get_smallest_values_of_interest(problem);
  }
  double upper() {
    return *coco_problem_get_largest_values_of_interest(problem);
  }
};
template <std::uint64_t Alg, int Dim, int Pop_Size, int Max> void coco_test() {
  using Alg_Type = sevobench::single_algorithm<Alg, Dim, Pop_Size, Max, false>;
  coco_suite_t *suite;
  coco_observer_t *observer;
  char suite_buffer[64];
  std::snprintf(suite_buffer, sizeof(suite_buffer), "dimensions:%d", Dim);
  suite = coco_suite("bbob", "", suite_buffer);
  char observe_buffer[128];
  std::snprintf(observe_buffer, sizeof(observe_buffer),
                "result_folder:Test,algorithm:%s", Alg_Type::name);
  observer = coco_observer("bbob", observe_buffer);
  coco_problem_t *problem;
  while ((problem = coco_suite_get_next_problem(suite, observer)) != NULL) {
    problem_wrapper p(problem);
    Alg_Type()(p, p.lower(), p.upper());
  }
  coco_observer_free(observer);
  coco_suite_free(suite);
}

#endif /* problem_wrapper_h */
