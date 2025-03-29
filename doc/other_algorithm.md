
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