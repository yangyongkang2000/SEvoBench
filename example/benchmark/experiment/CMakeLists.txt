cmake_minimum_required(VERSION 3.15)
get_filename_component(Data_Dir
        "${CMAKE_CURRENT_LIST_DIR}/../../../include/SEvoBench/problem/cec/"
        ABSOLUTE
)

if(NOT EXISTS ${Data_Dir})
    message(FATAL_ERROR "Data directory not found: ${Data_Dir}")
endif()
include("../../function.cmake")

find_package(ioh 0.3.18 QUIET)

add_executable(sevobench_benchmark_parallel sevobench_benchmark.cpp)
target_compile_definitions(sevobench_benchmark_parallel PRIVATE DATA_DIR="${Data_Dir}" PARALLEL=1)

add_executable(sevobench_benchmark sevobench_benchmark.cpp)
target_compile_definitions(sevobench_benchmark  PRIVATE DATA_DIR="${Data_Dir}" PARALLEL=0)

common(sevobench_benchmark_parallel)
common(sevobench_benchmark)

if(ioh_FOUND)
    add_executable(ioh_benchmark ioh_benchmark.cpp)
    target_link_libraries(ioh_benchmark PRIVATE ioh::ioh)
    common(ioh_benchmark)
endif()