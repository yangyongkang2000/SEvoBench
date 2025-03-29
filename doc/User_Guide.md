# SEvoBench User Guide
- [Introduction](#introduction)
    - [About SEvoBench](#about-sevobench)
    - [Installation](#installation)
        - [System Requirements](#system-requirements)
        - [Building the Project](#building-the-project)
        - [CMake](#using-cmake-for-package-management-in-sevobench)
## Introduction
### [About SEvoBench](#sevobench-user-guide)
&emsp;SEvoBench is a single-objective optimization algorithm testing framework written in C++. Its purpose is to facilitate single-objective evolutionary algorithm researchers to quickly and effectively test the optimization performance of algorithms. For this reason, the framework not only provides users with a series of interfaces for optimization algorithms and optimization problems, but also incorporates a certain number of state-of-the-art optimization algorithms and test sets.

### [Installation](#sevobench-user-guide)
&emsp;Simply download the source code from GitHub to your personal working directory. The computational kernel of SEvoBench is located in the *SEvoBench/inlcude/SEvoBench* folder, which contains all the necessary header files.

#### [System Requirements](#sevobench-user-guide)
    Compiler: Support for compilers with C++20 or above

#### [Building the Project](#sevobench-user-guide)
&emsp;As SEvoBench is implemented using C++ template programming, it does not need to be pre-compiled into a library. Simply include the project's header files and specify the project's header file search path.

#### [Using CMake for Package Management in SEvoBench](#sevobench-user-guide)
&emsp;SEvoBench supports package management using CMake.

&emsp;`cmake_minimum_required(VERSION 3.15)`

```shell
cmake -B build -S . -DBUILD_TESTS=OFF -DBUILD_EXAMPLE=OFF -DCMAKE_INSTALL_PREFIX=/* user-specified directory */   
cd build && cmake --build . && cmake --install .
```

+ [Compilation Options](#sevobench-user-guide)

&emsp;For GCC and Clang compilers, it is recommended to use the following compilation command:
```shell
-O3 -ffast-math -march=native -std=c++20 -lpthread
```
&emsp;For MSVC compiler, it is recommended to use the following compilation command:
```shell
/std:c++20 /O2 /fp:fast /arch:AVX2 /F 8388608
```
+ [Using CMake](#sevobench-user-guide)

```shell
cmake_minimum_required(VERSION 3.15)
project(quickstart LANGUAGES CXX)
set(EXE quickstart)
add_executable(${EXE} quickstart.cpp)
find_package(SEvoBench REQUIRED)
set_target_properties(${EXE} PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON)
target_link_libraries(${EXE} PUBLIC SEvoBench::SEvoBench)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CONFIGURATION_TYPES Release)
target_compile_options(${EXE} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /arch:AVX2 /F 8388608>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-ffast-math -march=native -Wall -Wextra>
)
target_link_libraries(${EXE} PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:pthread>
)
```
## Note on Version Migration and Namespace Changes

SEvoBench has undergone significant architectural revisions across its two major versions. In the first version, all built-in algorithm functions were organized under the namespace sevobench::other_algorithm. For users transitioning from the legacy version, please note the following:

### Backward Compatibility Reference:
The original implementation details remain accessible in the first version's documentation [other_algorithm.md](./other_algorithm.md).
This document serves as a historical reference but should be used with awareness of the updated namespace structure.
### Action Required for Migration:
Code referencing these algorithms must be updated to reflect the new namespace conventions in Version 2.
The revised architecture aims to improve modularity and maintainability, though it introduces breaking changes.
### Recommendation for New Development:
New implementations should adhere to the current version's namespace design.
Cross-check with the latest API documentation to ensure compatibility.

## Documentation Generation with Large Language Models**

The first version of SEvoBench required manual documentation writing, which imposed a significant workload for individual developers.
To address this challenge, the second major version has adopted large language models (LLMs) for comprehensive documentation generation.

The "./llm_code_input" directory contains all module codes provided to LLMs as input for documentation generation. These files are also made publicly available to serve as references for those who wish to use LLMs with SEvoBench data. Specifically:

1. [cec2017.md](./cec2017.md) and [problem.md](./problem.md) were generated using [llm_cec2017.txt](./llm_code_input/llm_cec2017.txt)
2. [de_module.md](./de_module.md) was produced from [llm_de.txt](./llm_code_input/llm_de.txt)
3. [pso_module.md](./pso_module.md) was created based on [llm_pso.txt](./llm_code_input/llm_pso.txt)
4. [experiment.md](./experiment.md) was generated using [llm_experiment.txt](./llm_code_input/llm_experiment.txt)

Important Notes:
- The LLM-generated documentation is provided for reference only
- We cannot guarantee the absolute accuracy of these documents
- Users are encouraged to refer to the original code in "./doc/llm_code_input" when using LLMs to generate their own documentation
- The provided input files may serve as useful templates for similar documentation tasks

This approach significantly reduces documentation workload while maintaining reasonable quality, though manual verification remains recommended for critical applications.

