# SEvoBench User Guide
- [Introduction](#introduction)  
     - [About SEvoBench](#about-sevobench)
     - [Installation](#installation)
       - [System Requirements](#system-requirements)
       - [Building the Project](#building-the-project)
       - [CMake](#using-cmake-for-package-management-in-sevobench)
    - [Getting Started](#getting-started)


- [Basics](#fundamentals)  
    - [Project Structure](#project-structure)
    - [Optimization Algorithm](#optimization-algorithms)
    - [Main Template Class of Algorithm](#main-template-class-of-algorithm)
    - [Test Problem Set](#test-problem-set)
    - [Main Template Class for Problem Sets](#main-template-class-for-problem-sets)
    - [Algorithm Test Functions](#algorithm-testing-function)
      - [About benchmark_result](#about-benchmark_resulttypename-t-typename-algorithm-typename-problem)
    - [Concepts and Constraints](#concepts-and-constraints)
- [API Reference](#api-introduction)
    - [APIs for Problem Sets](#about-the-problem-set-api)
      - [ShiftFunc](#shiftfunc)
      - [RotateShiftFunc](#rotateshiftfunc)
      - [AsyShiftFunc](#asyshiftfunc)
      - [AsyRotateShiftFunc](#asyrotateshiftfunc)
      - [CEC2020](#cec2020)
      - [CEC2022](#cec2022)
    - [APIs for Optimization Algorithms](#api-for-algorithms)
        - [Differential Evolution Algorithm](#differential-evolution-algorithm)
        - [Adaptive Differential Evolution Algorithm](#adaptive-differential-evolution-algorithm)
        - [Particle Swarm Optimization Algorithm](#particle-swarm-optimization-algorithm)
        - [Variant of Particle Swarm Optimization Algorithm for Large-Scale Optimization](#variants-of-particle-swarm-optimization-algorithms-for-solving-large-scale-optimization-problems)
        - [Biological-Inspired Optimization Algorithm](#bio-inspired-algorithms)
        - [Local Search Optimization Algorithm](#local-search-optimization-algorithms)
- [Examples](#examples)
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
### [Getting Started](#sevobench-user-guide)
&emsp;To demonstrate, we will use the built-in [DE algorithm](https://link.springer.com/article/10.1023/A:1008202821328) and the [CEC2020 test suite](https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm). For more detailed content, it will be explained later in the document.
+ The following code tests the [DE algorithm](https://link.springer.com/article/10.1023/A:1008202821328) on the [CEC2020 test suite](https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm).

```C++
#include"SEvoBench/parallel_algorithm_benchmark.hpp"
#include<cstdio>
int main() {
    auto data=sevobench::evo_bench<sevobench::DE, sevobench::CEC2020>();
    auto table=data.get_table_data();
    auto alg_name=data.algorithm_name();
    auto pro_name=data.problem_name();
    std::printf("Algorithm:%s Suite:%s,Dim:%d\n",alg_name,pro_name,data.problem_dim());
    for(int i=0;i<data.problem_size();i++) {
        std::printf("F%d mean:%f,std:%f,best:%f,time:%fms,1/4best:%f,middle:%f,3/4best:%f,worst:%f\n",i+1,table[8*i],table[8*i+1],table[8*i+2],table[8*i+3],table[8*i+4],table[8*i+5],table[8*i+6],table[8*i+7]);
    }
    return 0;
}
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

## [Fundamentals](#sevobench-user-guide)
### [Project Structure](#sevobench-user-guide)

| Content | Header File |
|---------|-------------|
| Base class for algorithms | single_algorithm.hpp |
| Base class for problems | single_problem.hpp |
| Parallel computation pool | parallel_task.hpp |
| Platform utility functions | tool.hpp |
| Optimization algorithms | pso.hpp, de.hpp, ... |
| Optimization problems | cec2020.hpp, cec2022.hpp, ... |
| Benchmark for optimization algorithms | parallel_algorithm_benchmark.hpp |
| Concept | common_concept.hpp|

### [Optimization Algorithms](#sevobench-user-guide)
&emsp;Since SEvoBench uses template programming, the optimization algorithms in the platform are implemented as template functions. If you are not familiar with template functions, please refer to [cppreference website](https://en.cppreference.com/w/cpp/language/function_template), and the book "C++ Templates: The Complete Guide" by David Vandevoorde and Nicolai M. Josuttis, or Chapter 16 of "C++ Primer, Fifth Edition".C++20 introduces [constraints and concepts](https://en.cppreference.com/w/cpp/language/constraints), which are also used in this framework to facilitate template programming.

&emsp;The typical declaration of an optimization algorithm in SEvoBench is as follows:
```C++
template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=algorithm_parameter<T>>
inline auto algorithm_func(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

```

&emsp;If you are not familiar with the *auto* keyword, please refer to [this link](https://en.cppreference.com/w/cpp/keyword/auto).

Now let me explain the meanings of the template parameters for the function:

+ *int Dim*: Integer type. It represents the dimension of the optimization problem function.

+ *int Pop_Size*: Integer type. For swarm intelligence-based optimization algorithms, it represents the initial population size of the algorithm. Users typically set a default value for a specific algorithm.

+ *int Max*: It represents the termination condition of the algorithm. When the template parameter *Memory_Flag* is True, it represents the maximum number of iterations for the algorithm. When *Memory_Flag* is False, it represents the maximum number of fitness evaluations for the optimization problem. Users typically set a default value for a specific algorithm.

+ *bool Memory_Flag*: Bool type. When *Memory_Flag* is True, the algorithm will record the convergence process and store the data in an array. When *Memory_Flag* is False, it will not record the convergence process. Users typically set a default value for a specific algorithm.

+ *typename F*: Specifies the optimization problem type. It does not need to be declared by users; it will be automatically deduced based on the function argument *function*.

+ *std::floating_point T*: Specifies the floating-point type required for calculations. It does not need to be declared by users; it will be automatically deduced based on the function argument.

+ *typename Parameter_Type*: Specifies the class that encapsulates the algorithm parameters. Algorithm parameters are crucial for an algorithm, but different algorithms have their own parameters. In order to have a unified interface for setting parameters, users need to pre-define a class to encapsulate the algorithm parameters and provide default values for them. For example, *algorithm_parameter\<T\>* can be used as the default type for *Parameter_Type*. Of course, *Parameter_Type* will be automatically deduced based on the function argument *pt*.

Now, let me explain the function arguments:

+ *F&& function*: For the function argument *F*, it should have the following characteristics:
```C++
auto result = function(x);
```
The input *x* of the *function* should be an array pointer, and *result* is the return value of the *function*.

+ *T l, T r*: SEvoBench assumes that the problem is a bound-constrained optimization problem, so *l* and *r* specify the lower and upper bounds of the variable domain.

+ *const Parameter_Type& pt = Parameter_Type()*: Users can use built-in classes or customize their own parameter classes, but make sure that the member variables of the class are consistent with the built-in classes.

&emsp;Due to the features of Modern C++, we can return different function value types based on different template parameters. When a user implements an *algorithm_func* function, SEvoBench requires the user to return different function values based on the template parameter type *Memory_Flag*.

```C++
template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=algorithm_parameter<T>>
inline auto algorithm_func(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept {
// middle process omitted
if constexpr(Memory_Flag)
        return std::make_tuple(best_pos, best_fit,convergence_curve);
if constexpr(!Memory_Flag)
    return std::make_pair(best_pos, best_fit);
}
```

+ SEvoBench assumes that the algorithm is seeking the global minimum of the problem.
+ When *Memory_Flag* is True, the function returns a [std::tuple](https://en.cppreference.com/w/cpp/utility/tuple) type that stores the optimal solution, optimal value, and array of convergence curve data.
+ When *Memory_Flag* is False, the function returns a [std::pair](https://en.cppreference.com/w/cpp/utility/pair) type that stores the optimal solution and optimal value of the problem.

&emsp;For the *convergence_curve* variable that records the convergence process, its data structure type is required to support the adaptation of the [std::begin()](https://en.cppreference.com/w/cpp/iterator/begin) and [std::end()](https://en.cppreference.com/w/cpp/iterator/end) functions. The storage format of the *convergence_curve* variable for storing convergence data is required as follows.

```C++
template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=algorithm_parameter<T>>
inline auto algorithm_func(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept {

/*
 Middle process omitted
*/

 /*
 Max_Itreation depends on the template variable Max, whether Max is the maximum number of iterations or the maximum number of fitness evaluations depends on the template variable B.
 */
for (int iteration=0; iteration < Max_Itreation; iteration++) {
    /*
    Middle process omitted
    */
    if constexpr (Memory_Flag) {
            convergence_curve[2*iteration] =Current_FES;// record the number of fitness evaluations after the current iteration
            convergence_curve[2*iteration+1] = Current_Best_Fit;// record the optimal fitness value after the current iteration
        }
 }

if constexpr(Memory_Flag)
        return std::make_tuple(best_pos, best_fit,convergence_curve);
if constexpr(!Memory_Flag)
    return std::make_pair(best_pos, best_fit);
}
```

&emsp;For each iteration, store the current number of fitness evaluations and the current best fitness value after each iteration in the *convergence_curve* variable.

### [Main Template Class of Algorithm](#sevobench-user-guide)

&emsp;In order to represent different algorithms with a unified algorithm class, SEvoBench uses the [partial specialization technique](https://en.cppreference.com/w/cpp/language/partial_specialization) to define a main template class for different optimization algorithm functions.


&emsp;For different optimization algorithm functions, a main template class is used to represent them. In the *single_algorithm.hpp* header file, namespace *sevobench* defines *single_algorithm* as the main template class.

```C++
namespace sevobench {

template<std::uint64_t Alg_HashName,int Dim,int Pop_Size,int Max,bool Memory_Flag>
struct single_algorithm ;

}
```

&emsp;Similarly, a main template class is used to represent different optimization algorithm parameters. In the *single_algorithm.hpp* header file, namespace *sevobench* defines *single_algorithm_parameter* as the main template class.

```C++
namespace sevobench {

template<std::uint64_t Alg_HashName,std::floating_point T>
struct single_algorithm_parameter;

}
```

+ The first template parameter *std::uint64_t Alg_HashName* represents the *Hash* value of the algorithm function. Therefore, each algorithm function has a unique *Alg_HashName* value, and the *Alg_HashName* value of different algorithm functions is different, so that the *single_algorithm* class can be partially specialized to express different algorithms with the main template *single_algorithm* class. The purpose of this design is to facilitate testing different algorithms on the test set.

+ The meanings of the other template parameters are the same as those of the template function parameters of the optimization algorithm.

&emsp;In the *single_algorithm.hpp* header file, there is also a macro expression defined called *REGISTER_ALGORITHM(func_name, func, parameter)*.

&emsp;The *REGISTER_ALGORITHM(func_name, func, parameter)* macro expression allows the algorithm function *func* to become a partial specialization of the *single_algorithm* class and generates a *func_name* variable representing the *Hash* value of the algorithm function *func*. After partial specialization, the *single_algorithm* class has a static member string array variable *name* (used to represent the algorithm name as *func_name*) and an [operator overload function](https://en.cppreference.com/w/cpp/language/operators) *inline auto operator()(F &&f, T l, T r)*. The specialized *single_algorithm* class can be treated as a [function object](https://en.cppreference.com/w/cpp/utility/functional). Similarly, the *parameter* generates a partial specialization of the *single_algorithm_parameter* class, which inherits from the parameter class. For more details, please refer to the implementation code in *single_algorithm.hpp* regarding the macro expression.

&emsp;Here is an example demonstrating how users can implement an algorithm and register it:

&emsp;First, include the necessary header files.
```C++
#include "single_algorithm.hpp"
/*
Include other header files required for algorithm implementation
*/
 
/* 
User-defined algorithm function Algorithm
*/
template<std::floating_point T>
struct Algorithm_Parameter {
    /*Definition of parameters required for the algorithm*/
};

template<int Dim, int Pop_Size, int Max, bool Memory_Flag, typename F, std::floating_point T, typename Parameter_Type = Algorithm_Parameter<T>>
inline auto Algorithm(F&& function, T l, T r, const Parameter_Type& pt = Parameter_Type()) noexcept {

/*
Declaration process omitted
*/

/*
Max_Itreation depends on template variable Max. Whether Max is the maximum number of iterations or the maximum number of fitness evaluations depends on the template variable B.
*/
    for (int iteration = 0; iteration < Max_Itreation; iteration++) {
        /*
        Omitted intermediate processes
        */
        if constexpr (Memory_Flag) {
            convergence_curve[2 * iteration] = Current_FES; // Record the number of fitness evaluations after the current iteration
            convergence_curve[2 * iteration + 1] = Current_Best_Fit; // Record the best fitness value after the current iteration
        }
    }

    if constexpr (Memory_Flag)
        return std::make_tuple(best_pos, best_fit, convergence_curve);
    if constexpr (!Memory_Flag)
        return std::make_pair(best_pos, best_fit);
}
```
Finally, the function is registered.

```C++
REGISTER_ALGORITHM(Algorithm_Name, Algorithm, Algorithm_Paramater)
```

In this way, users can implement the *Algorithm* function to obtain a partial specialization of the *single_algorithm* class and *single_algorithm_parameter* class. The partially specialized *single_algorithm_parameter* class inherits from *Algorithm_Paramater*, and also obtains a variable *Algorithm_Name* to represent the hash value of the *Algorithm* function.

We can use the partially specialized class to optimize the [Rastrigin problem](https://www.sfu.ca/~ssurjano/rastr.html):

```C++
#define ALGORITHM(ALG, Dim, F, L, R) single_algorithm<ALG, Dim, 30, 1000 * Dim, false>()(F, L, R)

int main() {
    std::cout << ALGORITHM(Algorithm_Name, 30, rastrigin<30>, -100.0f, 100.0f).second << "\n";
}
```

### [Test Problem Set](#sevobench-user-guide)

To perform continuous testing on a series of problems using an algorithm, it is necessary to construct a class for the problem set.

In SEvoBench, a problem set is generally declared using the following form:

```C++
template<int Prob_Index, int Dim, std::floating_point T>
struct Problem;
```

The design of the `Problem` class still utilizes [template partial specialization](https://en.cppreference.com/w/cpp/language/partial_specialization).

- Template parameter `int Prob_Index`: Represents the problem index.
- Template parameter `int Dim`: Represents the problem dimension.
- `std::floating_point T`: Floating-point type.

The first template parameter `int Prob_Index` is used to partially specialize the `Problem` class.

```C++
template<int Dim, std::floating_point T>
struct Problem<1, Dim, T> {
    static constexpr T L = -100; // Lower bound of the problem, using variable L
    static constexpr T U = 100; // Upper bound of the problem, using variable U
    
    /*
    Use operator overload function T operator()(T *x)
    */
    
    T operator()(T *x) noexcept {
        return 10 * Dim + std::accumulate(x, x + Dim, T(0), [](auto l, auto r) {
            return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
        });
    }
};
```
&emsp;We specialize a *Problem* class with an index of 1 and a problem type of [Rastrigin problem](https://www.sfu.ca/~ssurjano/rastr.html).

### [Main Template Class for Problem Sets](#sevobench-user-guide)
&emsp;To represent different problem sets with a unified problem set class, SEvoBench utilizes [template partial specialization](https://en.cppreference.com/w/cpp/language/partial_specialization).

&emsp;For different problem test sets, a main template class is used for representation.

&emsp;In the *single_problem.hpp* header file, the *sevobench* namespace declares *single_problem* as the main template class.

```C++
namespace sevobench {

template <std::uint64_t Prob_HashName, int Prob_Index, int Dim,
          std::floating_point T>
struct single_problem;

}
```

&emsp;For different problem sets, there are different initialization methods. Therefore, a main template class is also declared for initializing the problem sets.

```C++
namespace sevobench {

template <std::uint64_t Prob_HashName, int Dim, std::floating_point T>
struct init_single_problem;

}
```



+ Similar to the algorithm base class, the first template parameter of the class, *std::uint64_t Prob_HashName*, is the hash value of the problem set and is used for partial specialization.

+ The meaning of other template parameters is consistent with that of the test problem set.

&emsp;For this purpose, two macro expressions are provided:

`REGISTER_PROBLEM(problem_name, problem, size)`

`INIT_PROBLEM(problem, F)`

Used for partial specialization of problem sets and initialization methods.

&emsp;The following example shows how a user can implement a problem set and partially specialize it.

&emsp;First, include the necessary header files

```C++
#include"single_problem.hpp"


/*
the remaining header files needed for algorithm implementation
*/
```

&emsp;Declare and define a problem, define an initialization method.
```C++
template< int Prob_Index, int Dim,std::floating_point T >
Struct  Problem;



template<int Dim, std::floating_point T>
struct Problem<1, Dim, T> {
    static constexpr T L = -100; // Lower bound of the problem, using variable L
    static constexpr T U = 100; // Upper bound of the problem, using variable U
    
    /*
    Use operator overload function T operator()(T *x)
    */
    
    T operator()(T *x) noexcept {
        return 10 * Dim + std::accumulate(x, x + Dim, T(0), [](auto l, auto r) {
            return l + r * r - 10 * std::cos(2 * std::numbers::pi_v<T> * r);
        });
    }
};

template<std::floating_point T,int Dim>
void init_Problem() {
    /*
    user-defined
    */
}
```
&emsp;Partial specialization of problem and initialization method



```C++
REGISTER_PROBLEM(Problem_Name, Problem, 1) // The last parameter represents the number of problems in the problem set
INIT_PROBLEM(Problem, init_Problem)
```



&emsp;By implementing the `Problem` problem set and partially specializing it as a `single_problem` class, users can obtain a variable `Problem_Name` to represent the hash value of the `Problem` problem set.

&emsp;However, in general, users do not need to add problem sets themselves, as SEvoBench platform already includes a series of built-in problem sets.

### [Algorithm Testing Function](#sevobench-user-guide)
In the header file *parallel_algorithm_benchmark.hpp*, under the namespace *sevobench*, we have defined the following template function:

```C++
namespace sevobench {

template<std::uint64_t Alg_HashName,std::uint64_t Prob_HashName,bool Memory_Flag,int Dim,int Pop_Size,int Max,int Runs,bool parallel,std::floating_point T,typename Parameter_Type=typename single_algorithm_parameter<Alg_HashName, T>::base_type>
inline auto evo_bench(const Parameter_Type & par=Parameter_Type()) noexcept;

}
```

For users, it is only necessary to pay attention to the following template parameters:

+ *std::uint64_t Alg_HashName*: The hash value specifying a specific algorithm.
+ *std::uint64_t Prob_HashName*: The hash value specifying a specific problem set.
+ *bool Memory_Flag*: Whether to record the convergence process. The default value is false.
+ *int Dim*: The dimension of the problem set. The default value is 30.
+ *int Pop_Size*: The population size for the evolutionary algorithm. The default value is 30.
+ *int Max*: The maximum number of fitness evaluations or algorithm iterations for the problem. The default value is 1000 times the dimension.
+ *int Runs*: The number of test repetitions. The default value is 30.
+ *bool parallel*: Whether to perform parallel computation. The default value is true.
+ *std::floating_point T*: The floating-point type. The default value is float.
+ *typename Parameter_Type*: Specifies the class that encapsulates the algorithm parameters. The default type is Algorithm_Paramater from the macro *REGISTER_ALGORITHM(Algorithm_Name,Algorithm,Algorithm_Paramater)*.

The testing function returns a value of type *benchmark_result<std::floating_point T,typename Algorithm,typename Problem>*. This data structure stores the results of the algorithm testing.

#### [About *benchmark_result\<std::floating_point T, typename Algorithm, typename Problem>*](#sevobench-user-guide)

+ Defined in the header file *parallel_algorithm_benchmark.hpp*, under the namespace *sevobench::sevobench_detail*.

+ Declared as follows:

```C++
namespace sevobench::sevobench_detail {

template<std::floating_point T,typename Algorithm,typename Problem>
class benchmark_result;

}
```

+ The member functions of this class include:

> *auto get_graph_data()* : Returns the convergence data *graph_data*, stored in one-dimensional array form with data structure *std::vector\<T\>*.

> *auto get_table_data()* : Returns the statistical data *table_data*, stored in one-dimensional array form with data structure *std::vector\<T\>*.

> *constexpr auto problem_dim()* : Returns the dimension of the problem set being tested.

> *constexpr auto problem_size()* : Returns the number of problems in the problem set being tested.

> *constexpr auto problem_name()* : Returns the name of the problem set being tested.

> *constexpr auto algorithm_name()* : Returns the name of the algorithm being tested.

> *constexpr auto algorithm_max()* : Returns the maximum number of iterations for the algorithm or the maximum number of fitness evaluations for the problem.

> *constexpr auto algorithm_flag()* : Returns the value of the *Memory_Flag* flag for the algorithm.

 > *void write_table(const char \*dir_name)* : Writes the statistical data to a file named `{Algorithm}_{Populations}_{Problem}_{Dim}_table.txt` in the directory *dir_name*.

Each line of the file represents a piece of data, and every 8 pieces of data correspond to one group of results for a specific problem in the problem set. These 8 pieces of data represent the mean value, standard deviation, best value, average time, lower quartile, median, upper quartile and worst value for the algorithm's performance on the problem.

> *void write_curve(const char \*dir_name)* : Writes the convergence data to a file named `{Algorithm}_{Populations}_{Problem}_F{index}_{Dim}_curve.txt` in the directory *dir_name*. Multiple files will be generated, each corresponding to a specific problem in the problem set and a specific run of the algorithm.

Each line of the file represents a point in the convergence curve of the algorithm for the corresponding problem. The point consists of the number of fitness evaluations performed (*FES*) and the best fitness value obtained so far (*Best_Fit*).

+ It should be noted that the data types returned by the functions *auto get_graph_data()* and *auto get_table_data()* have the same layout as those in the text files.

## [Concepts and Constraints](#sevobench-user-guide)
Some common concepts are defined in *common_concept.hpp* to constrain template parameters in algorithm function declarations.

```C++
namespace sevobench {

template <auto Dim, auto Pop_Size, auto Max>
concept algorithm_parameter_concept 

template <typename F, typename T>
concept algorithm_problem_concept 

template <typename F, typename G, typename T>
concept algorithm_positions_concept 

template <typename F, typename G, typename T>
concept algorithm_vector_concept 

template <auto Dim, auto Pop_Size, auto Max, typename F, typename T>
concept algorithm_func_concept

}
```



## [API Introduction](#sevobench-user-guide)

### [About the Problem Set API](#sevobench-user-guide)

#### [*ShiftFunc*](#sevobench-user-guide)

+ Hash variable name of the problem set: *sevobench::ShiftFunc*
+ Problem set type: Classic single-objective functions with 13 expandable dimensions, each function undergoes random shift processing
+ Related literature link for the problem set: [https://ieeexplore.ieee.org/document/771163](https://ieeexplore.ieee.org/document/771163)
+ Located in header file: *classic_problem.hpp single_problem.hpp*

&emsp; Relevant declaration under *namespace sevobench* in *single_problem.hpp*
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<ShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*RotateShiftFunc*](#sevobench-user-guide)

+ Hash variable name of the problem set: *sevobench::RotateShiftFunc*
+ Problem set type: Classic single-objective functions with 13 expandable dimensions, each function undergoes random rotation shift processing
+ Related literature link for the problem set: [https://ieeexplore.ieee.org/document/771163](https://ieeexplore.ieee.org/document/771163)
+ Located in header file: *classic_problem.hpp single_problem.hpp*

&emsp; Relevant declaration under *namespace sevobench* in *single_problem.hpp*
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<RotateShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*AsyShiftFunc*](#sevobench-user-guide)

+ Hash variable name of the problem set: *sevobench::AsyShiftFunc*
+ Problem set type: Classic single-objective functions with 13 expandable dimensions, each function undergoes random shift processing and non-linear transformations to exhibit ill-conditioning, irregularity, and symmetry-breaking characteristics.
+ Related literature link for the problem set: [https://ieeexplore.ieee.org/document/771163](https://ieeexplore.ieee.org/document/771163)
+ Related literature link for the function transformations: Dim. Hansen, S. Finck, R. Ros, A. Auger, Real-Parameter Black-Box Optimization Benchmarking 2009: Noiseless Functions Definitions, Tech. rep. RR-6829, INRIA, 2010.
+ Located in header file: *classic_problem.hpp single_problem.hpp*

&emsp; Relevant declaration under *namespace sevobench* in *single_problem.hpp*
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<AsyShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*AsyRotateShiftFunc*](#sevobench-user-guide)

+ Hash variable name of the problem set: *sevobench::AsyRotateShiftFunc*
+ Problem set type: Classic single-objective functions with 13 expandable dimensions, each function undergoes random rotation shift processing and non-linear transformations to exhibit ill-conditioning, irregularity, and symmetry-breaking characteristics.
+ Related literature link for the problem set: [https://ieeexplore.ieee.org/document/771163](https://ieeexplore.ieee.org/document/771163)
+ Related literature link for the function transformations: Dim. Hansen, S. Finck, R. Ros, A. Auger, Real-Parameter Black-Box Optimization Benchmarking 2009: Noiseless Functions Definitions, Tech. rep. RR-6829, INRIA, 2010.
+ Located in header file: *classic_problem.hpp single_problem.hpp*

&emsp; Relevant declaration under *namespace sevobench* in *single_problem.hpp*
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<AsyRotateShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*CEC2020*](#sevobench-user-guide)

+ Hash variable name of the problem set: *sevobench::CEC2020*
+ Problem set type: 10 single-objective functions from IEEE CEC 2020 benchmark, with expandable dimensions. The original benchmark limits the dimensions of the functions, while the benchmark in this platform supports expandable dimensions.
+ Related literature link for the problem set: [https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm](https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm)
+ Located in header file: *cec2020.hpp single_problem.hpp*

&emsp; Relevant declaration under *namespace sevobench* in *single_problem.hpp*
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<CEC2020,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*CEC2022*](#sevobench-user-guide)

+ Hash variable name of the problem set: *sevobench::CEC2022*
+ Problem set type: 12 single-objective functions from IEEE CEC 2022 benchmark, with expandable dimensions. The original benchmark limits the dimensions of the functions, while the benchmark in this platform supports expandable dimensions.
+ Related literature link for the problem set: [https://www3.ntu.edu.sg/home/epnsugan/index_files/CEC2022/CEC2022.htm](https://www3.ntu.edu.sg/home/epnsugan/index_files/CEC2022/CEC2022.htm)
+ Located in header file: *cec2022.hpp single_problem.hpp*

&emsp; Relevant declaration under *namespace sevobench* in *single_problem.hpp*
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<CEC2022,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### Note

&emsp;Due to the random generation of relevant test data during each problem set initialization, it is recommended to compare the performance of different algorithms on the same problem set initialization.

### [API for Algorithms](#sevobench-user-guide)

#### [Differential Evolution Algorithm](#sevobench-user-guide)
+ Hash variable name of the algorithm: *sevobench::DE*
+ Related literature link for the algorithm: https://link.springer.com/article/10.1023/A:1008202821328
+ Header file location: *de.hpp*
+ Function declaration:
```C++
namespace sevobench {

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=de_parameter<T>>
inline auto de_optimize(G &&positions,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=de_parameter<T>>
inline auto de(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```
+ Parameter class *de_parameter*: *f* represents the weight coefficient, and *cr* represents the crossover probability.
```C++
namespace sevobench {

template<std::floating_point T>
struct de_parameter {
    T f=0.5;
    T cr=0.9;
};

}
```

<br>
<br>

#### [Adaptive Differential Evolution Algorithm](#sevobench-user-guide)
+ Hash variable name of the algorithm: *sevobench::JADE*
+ Related literature link for the algorithm: https://ieeexplore.ieee.org/document/5208221
+ Header file location: *jade.hpp*
+ Function declaration:
```C++
namespace sevobench {

    template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=jade_parameter<T>>
inline auto jade_optimize(G &&positions,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type());

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=jade_parameter<T>>
inline auto jade(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;
 
}
```
+ Parameter class *jade_parameter*: The parameters are consistent with the paper.
```C++
namespace sevobench {

template<std::floating_point T>
struct jade_parameter {
    T c=0.1;
    T p=0.1;
};

}
```
<br>
<br>

+ Hash variable name of the algorithm: *sevobench::SHADE*
+ Related literature link for the algorithm: https://ieeexplore.ieee.org/document/6557555  https://ieeexplore.ieee.org/abstract/document/6900380 
+ Header file location: *template_shade.hpp*
+ Function declaration:
```C++
namespace sevobench {

   template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=shade_parameter<T>>
inline auto shade(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept ;

}
```

<br>
<br>

+ Hash variable name of the algorithm: *sevobench::LSHADE*
+ Related literature link for the algorithm: https://ieeexplore.ieee.org/abstract/document/6900380 
+ Header file location: *template_shade.hpp*
+ Function declaration:
```C++
namespace sevobench {

template<int Dim,int Pop_Size=18*Dim,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=lshade_parameter<T>>
inline auto lshade(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ Parameter classes *shade_parameter* and *lshade_parameter*: The parameters are consistent with the paper.

```C++
namespace sevobench {

template<std::floating_point T>
struct shade_parameter {
    T r_Arc{2};
    T p=T(0.1);
    int H{100};
};

template<std::floating_point T>
struct lshade_parameter {
    T r_Arc{2.6};
    T p=T(0.11);
    int H{6};
};

}
```

<br>
<br>

#### [Particle Swarm Optimization Algorithm](#sevobench-user-guide)

+ Algorithm hash variable name: `sevobench::PSO`
+ Algorithm related literature link: [https://ieeexplore.ieee.org/document/699146](https://ieeexplore.ieee.org/document/699146)
+ Header file: `pso.hpp`
+ Function declarations:

```C++
namespace sevobench {

template<int Dim, int Pop_Size = 30, int Max = 1000*Dim, bool Memory_Flag = false, typename G, typename F, std::floating_point T, typename Parameter_Type=pso_parameter<T>>
inline auto pso_optimize(G &&positions, F &&f, T left_bound, T right_bound, const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim, int Pop_Size = 30, int Max = 1000*Dim, bool Memory_Flag = false, typename F, std::floating_point T, typename Parameter_Type=pso_parameter<T>>
inline auto pso(F &&f, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ Algorithm parameter `pso_parameter`:

```C++
namespace sevobench {

template<std::floating_point T>
struct pso_parameter {
    T c1 = T(2);
    T c2 = T(2);
    T w_min = T(0.4);
    T w_max = T(0.9);
    T v_max_ratio = T(0.2);
};

}
```

+ Algorithm hash variable name: `sevobench::SPSO2007`
+ Algorithm related literature link: [https://ieeexplore.ieee.org/document/4223164](https://ieeexplore.ieee.org/document/4223164), [https://hal.science/hal-00764996](https://hal.science/hal-00764996)
+ Header file: `template_pso.hpp`
+ Function declaration:

```C++
namespace {

template<int Dim, int Pop_Size=30, int Max=1000*Dim, bool Memory_Flag=false, typename F, std::floating_point T, typename Parameter_Type=spso2007_parameter<T>>
inline auto spso2007(F&& function, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ Algorithm hash variable name: `sevobench::SPSO2011`
+ Algorithm related literature link: [https://ieeexplore.ieee.org/document/6557848](https://ieeexplore.ieee.org/document/6557848), [https://hal.science/hal-00764996](https://hal.science/hal-00764996)
+ Header file: `template_pso.hpp`
+ Function declaration:

```C++
namespace sevobench {

template<int Dim, int Pop_Size=30, int Max=1000*Dim, bool Memory_Flag=false, typename F, std::floating_point T, typename Parameter_Type=spso2011_parameter<T>>
inline auto spso2011(F&& function, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ Algorithm parameters `spso2007_parameter` and `spso2011_parameter`: The parameters are consistent with the parameters in the literature [https://hal.science/hal-00764996](https://hal.science/hal-00764996).

```C++
namespace sevobench {

template<std::floating_point T>
struct spso2007_parameter {
    T w = T(1) / T(2 * std::numbers::ln2_v<T>);
    T c = T(0.5) + std::numbers::ln2_v<T>;
};

template<std::floating_point T>
struct spso2011_parameter {
    T w = T(1) / T(2 * std::numbers::ln2_v<T>);
    T c = T(0.5) + std::numbers::ln2_v<T>;
};

}
```


<br>
<br>




#### [Variants of Particle Swarm Optimization Algorithms for Solving Large-Scale Optimization Problems](#sevobench-user-guide)

+ Algorithm hash variable name: `sevobench::CSO`
+ Algorithm related literature link: [https://ieeexplore.ieee.org/document/6819057](https://ieeexplore.ieee.org/document/6819057)
+ Header file: `cso.hpp`
+ Function declarations:

```C++
namespace sevobench {

template<int Dim, int Pop_Size = 250, int Max = 1000*Dim, bool Memory_Flag = false, typename G, typename F, std::floating_point T, typename Parameter_Type=cso_parameter<T>>
inline auto cso_optimize(G &&positions, F&& function, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim, int Pop_Size = 250, int Max = 1000*Dim, bool Memory_Flag = false, typename F, std::floating_point T, typename Parameter_Type=cso_parameter<T>>
inline auto cso(F&& function, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ Algorithm parameter `cso_parameter`:The parameters are consistent with the paper.

```C++
namespace sevobench {

template<std::floating_point T>
struct cso_parameter {
    T fai = T(0.1);
};

}
```

+ Algorithm hash variable name: `sevobench::SLPSO`
+ Algorithm related literature link: [https://www.sciencedirect.com/science/article/pii/S0020025514008366](https://www.sciencedirect.com/science/article/pii/S0020025514008366)
+ Header file: `slpso.hpp`
+ Function declarations:

```C++
namespace sevobench {

template<int Dim, int Pop_Size, int Max, bool Memory_Flag = false, typename G, typename F, std::floating_point T, typename Parameter_Type=slpso_parameter<T>>
inline auto slpso_optimize(G &&positions, F&& function, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim, int Pop_Size, int Max, bool Memory_Flag = false, typename F, std::floating_point T, typename Parameter_Type=slpso_parameter<T>>
inline auto slpso(F&& function, T l, T r, const Parameter_Type& pt=Parameter_Type()) noexcept ;

}
```

+ Algorithm parameter `slpso_parameter`:The parameters are consistent with the paper.

```C++
namespace sevobench {

template<std::floating_point T>
struct slpso_parameter {
    T beta = T(0.01);
    T alpha = T(0.5);
};

}
```
<br>
<br>

#### [Bio-inspired Algorithms](#sevobench-user-guide)
+ Artificial Bee Colony Algorithm
+ Hash variable name: *sevobench::ABC*
+ Relevant literature:

> D. Karaboga, An Idea Based on Honey Bee Swarm for Numerical Optimization, Technical Report-TR06, Erciyes University, Engineering Faculty, Computer Engineering Department 2005.

> D. Karaboga, B. Basturk, A Powerful and Efficient Algorithm for Numerical Function Optimization: Artificial Bee Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39, Issue:3,pp:459-171, November 2007,ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x 

> D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony (ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, Pages 687-697. 

> D. Karaboga, B. Akay, A Comparative Study of Artificial Bee Colony Algorithm,  Applied Mathematics and Computation, 214, 108-132, 2009. 

+ Header file location: *abc.hpp*
+ Function declaration:
```C++
namespace sevobench {

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=abc_parameter<T>>
inline auto abc_optimize(G&& sol,F &&f,T d,T w,const Parameter_Type &pt=Parameter_Type()) noexcept ;

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=abc_parameter<T>>
inline auto abc(F &&f,T d,T w,const Parameter_Type &pt=Parameter_Type()) noexcept;

}
```

+ Algorithm parameter: *abc_parameter*
```C++
namespace sevobench {

template<std::floating_point T>
struct abc_parameter {
    T ratio=T(0.5);
};

}
```
<br>
<br>

#### [Local Search Optimization Algorithms](#sevobench-user-guide)

+ Random Search Algorithm
+ Hash variable name: *sevobench::RandomSearch*
+ Header file location: *random_search.hpp*
+ Function declaration:
```C++
namespace sevobench {

template<int Dim,int Pop_Size=1,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type>
inline auto random_search_optimize(G&& sol,F&& function,T l,T r,const Parameter_Type&) noexcept ;

template<int Dim,int Pop_Size=1,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type>
inline auto random_search(F&& function,T l,T r,const Parameter_Type&pt) noexcept ;

}
```

<br>
<br>

+ (1+1) Evolution Strategy
+ Hash variable name: *sevobench::ES*
+ Relevant literature:

> Rechenberg, I. (1973). Evolutionsstrategie: Optimierung technischer Systeme nach Prinzipien der biologischen Evolution. Stuttgart: Fromann-Holzboog.

> Beyer, H.-G., & Schwefel, H.-P. (2002). Evolution strategies - A comprehensive introduction. Natural computing, 1(1), 3-52.

+ Header file location: *es.hpp*
+ Function declaration:
```C++
namespace sevobench {

template<int Dim,int Pop_Size=1,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=es_parameter<T>>
inline auto es_optimize(G&& sol,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim,int Pop_Size=1,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=es_parameter<T>>
inline auto es(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ Algorithm parameter: *es_parameter*
```C++
namespace sevobench {

template<std::floating_point T>
struct es_parameter {
    T sigma=T(1);
};

}
```

<br>
<br>

## [Examples](#sevobench-user-guide)
For user convenience, five examples are provided in the *example* folder:
+ *example/quickstart*: Demonstrates how to quickly get started with the SEvoBench framework.
+ *example/coco_test*: Illustrates how to integrate SEvoBench algorithms with the BBOB test suite on the COCO platform.
+ *example/cmaes_test*: Demonstrates the integration of the open-source CMA-ES algorithm code from GitHub into the SEvoBench framework.
+ *example/parallel_test*: Showcases the parallel testing performance of SEvoBench.
+ *example/add_problem*: Illustrates how to add a problem to the SEvoBench framework.