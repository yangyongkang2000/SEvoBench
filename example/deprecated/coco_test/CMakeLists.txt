cmake_minimum_required(VERSION 3.15)
project(coco_test LANGUAGES CXX C)
include("../function.cmake")
add_library(coco STATIC coco.c)
add_executable(coco_test
        main.cpp coco.c)
common(coco_test)
common(coco)