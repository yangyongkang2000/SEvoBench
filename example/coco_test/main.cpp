//
//  main.cpp
//  coco_test
//
//  Created by 杨永康 on 2023/9/11.
//

#include "problem_wrapper.hpp"
#include <iostream>

int main() {
  // insert code here...
  coco_test<sevobench::SHADE, 40, 100, 400000>();
  coco_test<sevobench::JADE, 40, 100, 400000>();
  coco_test<sevobench::LSHADE, 40, 720, 400000>();
  coco_test<sevobench::SLPSO, 40, 100, 400000>();
  coco_test<sevobench::SPSO2011, 40, 40, 400000>();
  coco_test<sevobench::SPSO2007, 40, 25, 400000>();
  return 0;
}
