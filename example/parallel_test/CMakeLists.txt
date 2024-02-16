cmake_minimum_required(VERSION 3.15)
project(parallel_test)
include("../function.cmake")
add_executable(parallel
        main.cpp)
add_executable(single main.cpp)

common(parallel)

common(single)

target_compile_definitions(parallel PRIVATE PARALLEL=true)
target_compile_definitions(single PRIVATE PARALLEL=false)
