cmake_minimum_required(VERSION 3.15)
project(tests)
enable_testing()
include("../example/function.cmake")
add_executable(test_tool test_tool.cpp)
add_executable(test_problem test_problem.cpp)
add_executable(test_parallel_task test_parallel_task.cpp)
add_executable(test_algorithm test_algorithm.cpp)
add_executable(test_benchmark test_benchmark.cpp)
add_executable(test_positions test_positions.cpp)
common(test_tool)
common(test_problem)
common(test_parallel_task)
common(test_algorithm)
common(test_benchmark)
common(test_positions)
add_test(NAME test_tool COMMAND test_tool)
add_test(NAME test_problem COMMAND test_problem)
add_test(NAME test_algorithm COMMAND test_algorithm)
add_test(NAME test_benchmark COMMAND test_benchmark)
add_test(NAME test_parallel_task COMMAND test_parallel_task)
add_test(NAME test_positions COMMAND test_positions)