# SEvoBench用户指南

- [介绍](#introduction)  
     - [关于SEvoBench](#关于sevobench)
     - [安装](#安装)
       - [系统要求](#系统要求)
       - [编译项目](#编译项目)
       - [使用CMake进行包管理](#使用cmake进行包管理)
    - [快速开始](#快速开始)


- [基础](#basics)  
    - [项目结构](#项目结构)
    - [优化算法](#优化算法)
    - [算法主模版类](#算法主模版类)
    - [测试问题集](#测试问题集)
    - [问题集主模版类](#问题集主模版类)
    - [算法测试函数](#算法测试函数)
      - [关于benchmark_result](#关于benchmark_resulttypename-ttypename-algorithmtypename-problem)
    - [概念与约束](#概念与约束)
- [API介绍](#api介绍)
    - [关于问题集的API](#关于问题集的api)
      - [ShiftFunc](#shiftfunc)
      - [RotateShiftFunc](#rotateshiftfunc)
      - [AsyShiftFunc](#asyshiftfunc)
      - [AsyRotateShiftFunc](#asyrotateshiftfunc)
      - [CEC2020](#cec2020)
      - [CEC2022](#cec2022)
    - [关于算法的API](#关于算法的api)
        - [差分进化算法](#差分进化算法)
        - [自适应差分进化算法](#自适应差分进化算法)
        - [粒子群优化算法](#粒子群优化算法)
        - [用于解决大规模优化的粒子群优化算法变种](#用于解决大规模优化的粒子群优化算法变种)
        - [生物启发式算法](#生物启发式算法)
        - [局部搜索优化算法](#局部搜索优化算法)
- [示例](#示例)



## 介绍
### [关于SEvoBench](#sevobench用户指南)
&emsp;SEvoBench是一个使用C++编写的单目标优化算法测试框架.它的目的是方便单目标进化算法研究工作者快速有效测试算法的优化性能,为此框架不仅为用户提供了优化算法和优化问题的一系列接口,而且还内置了一定数量的state of the art 的优化算法和测试集.

### [安装](#sevobench用户指南)
&emsp;从GitHub上下载源码到用户个人工作文件夹下即可.SEvoBench计算内核位于*SEvoBench/include/SEvoBench*文件夹下,该文件夹包含了所需要的所有头文件.

#### [系统要求](#sevobench用户指南)
    Compile:支持C++20以上的编译器

#### [编译项目](#sevobench用户指南)
&emsp;由于SEvoBench基于C++模版编程实现,无需提前编译成链接库,只需引入项目头文件,指定项目头文件搜索路径即可.

#### [使用CMake进行包管理](#sevobench用户指南)
&emsp;SEvoBench支持使用CMake进行包管理

&emsp;cmake_minimum_required(VERSION 3.15)
```shell
cmake -B build  -S . -DBUILD_TESTS=OFF -DBUILD_EXAMPLE=OFF -DCMAKE_INSTALL_PREFIX=/* 用户指定文件夹 */   
cd build &&cmake --build .&& cmake --install .
```

### [快速开始](#sevobench用户指南)
&emsp;为了演示我们使用内置的[DE算法](https://link.springer.com/article/10.1023/A:1008202821328)和[CEC2020测试集](https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm).有关更详细的内容,在文档后面会讲解.
+ 以下代码为[DE算法](https://link.springer.com/article/10.1023/A:1008202821328)在[CEC2020测试集](https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm)上测试
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


+ [编译选项](#sevobench用户指南)

&emsp;对于GCC Clang编译器,推荐编译命令
```shell
-O3 -ffast-math -march=native -std=c++20 -lpthread
```
&emsp;对于MSVC编译器,推荐编译命令
```shell
/std:c++20 /O2 /fp:fast  /arch:AVX2 /F 8388608
```

+ [使用CMake](#sevobench用户指南)

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

## [基础](#sevobench用户指南)
### [项目结构](#sevobench用户指南)

| 内容 | 头文件|
|-----|------|
|  algorithm基类 | single_algorithm.hpp|
|  problem基类 | single_problem.hpp|
| 并行计算池 | parallel_task.hpp|
| 平台工具函数 | tool.hpp |
| 优化算法 | pso.hpp de.hpp ...|
| 优化问题 | cec2020.hpp cec2022.hpp ...|
| 优化算法测试| parallel_algorithm_benchmark.hpp|
| 约束与概念     | common_concept.hpp |


### [优化算法](#sevobench用户指南)
&emsp;由于SEvoBench使用模版编程的风格,所以平台中的优化算法使用模版函数风格,如果你对模版函数不了解,请参阅[cppreference网站](https://en.cppreference.com/w/cpp/language/function_template), David Vandevoorde and Nicolai M. Josuttis.所著的《 C++ template: the complete guide》.《C++ Primer Fifth Edition》第十六章内容.
C++20引入了[约束与概念](https://en.cppreference.com/w/cpp/language/constraints),本框架也使用[约束与概念](https://en.cppreference.com/w/cpp/language/constraints)用于更好的进行模版编程.

&emsp;SEvoBench中优化算法一般使用如下声明
```C++
template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=algorithm_parameter<T>>
inline auto algorithm_func(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

```
&emsp;如果你不了解*auto*关键词,[请参阅](https://en.cppreference.com/w/cpp/keyword/auto)


&emsp;解释一下函数的模版参数代表的意义

+ *int Dim*:    整数类型.指代优化问题函数的维度.

+ *int Pop_Size*:    整数类型.对于群智能优化算法来说,指代算法初始种群的大小.用户一般对特定算法设置默认值.

+ *int Max*:  指代算法的终止条件.当模版参数*Memory_Flag*为True时,指代算法最大迭代次数.当模版参数*Memory_Flag*为False时.指代优化问题的适应度最大计算次数.用户一般对特定算法设置默认值.

+ *bool Memory_Flag*:   Bool类型.当*Memory_Flag*为True时,算法计算过程中,会记录收敛过程,用数组来存储数据;当*Memory_Flag*为False时,不会记录收敛过程.用户一般对特定算法设置默认值.

+ *typename F*:指定优化问题类型.不需要用户去声明,会根据函数实参*function*自动推导.

+ *std::floating_point T*:指定计算所需的浮点类型.不需要用户去声明,会根据函数实参自动推导.

+ *typename Parameter_Type*:指定封装算法参数的类.算法参数作为算法很重要的一部分,但是不同算法有着各自的参数,为了对参数设置有统一的接口,用户需要提前将算法参数封装成一个类,并为参数设置默认值.比如将封装好的  *algorithm_parameter\<T\>*  作为 *Parameter_Type* 默认类型,当然*Parameter_Type*也会根据函数实参*pt*自动推导

对函数实参的说明

+ *F&& function*:对于函数实参的*F*类型*function*,要求具有以下特性
```C++
auto result=function(x);// function的实参x要求为数组类型指针,result为function函数计算返回的结果.
``` 

+ *T l,T r*:SEvoBench假设问题为Bound constrained optimization problems,所以*l*,*r*分别指定问题变量域的下界和上界.

+ *const Parameter_Type& pt=Parameter_Type()*:用户可以使用内置的类,也可自定义参数类,但确保类成员变量和内置的类保持一致.




&emsp;由于Modern C++的特性,我们可以根据模版参数的不同,返回不同的函数值类型.当用户自己实现一个*algorithm_func*函数时,SEvoBench要求用户根据模版参数类型*Memory_Flag*的不同返回不同的函数值.
```C++
template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=algorithm_parameter<T>>
inline auto algorithm_func(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept {
// 中间过程省略
if constexpr(Memory_Flag)
        return std::make_tuple(best_pos, best_fit,convergence_curve);
if constexpr(!Memory_Flag)
    return std::make_pair(best_pos, best_fit);
}
```
+ SEvoBench假定算法求解问题的全局最小值
+ 当*Memory_Flag*为True时,函数返回[std::tuple](https://en.cppreference.com/w/cpp/utility/tuple)类型,分别存储问题的最优解,最优值,记录收敛过程的数组.
+ 当*Memory_Flag*为False时,函数返回[std::pair](https://en.cppreference.com/w/cpp/utility/pair)类型,分别存储问题的最优解,最优值.

&emsp;对记录收敛过程的*convergence_curve*变量,要求其数据结构类型支持适配[std::begin()](https://en.cppreference.com/w/cpp/iterator/begin) 和[std::end()](https://en.cppreference.com/w/cpp/iterator/end)函数.对于*convergence_curve*变量如何存储收敛数据,要求如下形式存储.
```C++
 template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=algorithm_parameter<T>>
inline auto algorithm_func(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept {

/*
 声明过程省略
 */

  /*
  Max_Itreation取决于模版变量Max,Max是最大迭代次数还是适应度最大计算次数则取决于模版变量B
  */
 for (int iteration=0; iteration < Max_Itreation; iteration++) {
    /*
    中间过程省略
    */
    if constexpr (Memory_Flag) {
            convergence_curve[2*iteration] =Current_FES;//记录当前迭代之后的问题适应度计算次数
            convergence_curve[2*iteration+1] = Current_Best_Fit;//记录当前迭代之后的最佳适应度值
        }
 }


if constexpr(Memory_Flag)
        return std::make_tuple(best_pos, best_fit,convergence_curve);
if constexpr(!Memory_Flag)
    return std::make_pair(best_pos, best_fit);
}
```
&emsp;对于每一次迭代,存储每一次迭代之后的当前适应度计算次数和当前最佳适应度值.


### [算法主模版类](#sevobench用户指南)
&emsp;为了使得不同的算法有统一的算法类去表示,SEvoBench使用[模版偏特化技术](https://en.cppreference.com/w/cpp/language/partial_specialization) 

&emsp;对不同的优化算法函数,用一个主模版类去表示. 

&emsp;在*single_algorithm.hpp*头文件中namespace *sevobench*下定义了*single_algorithm*作为主模版类
```C++
namespace sevobench {

template<std::uint64_t Alg_HashName,int Dim,int Pop_Size,int Max,bool Memory_Flag>
struct single_algorithm ;

}
```


&emsp;同理,对不同的优化算法参数,也用一个主模版类去表示.在*single_algorithm.hpp*头文件中在namespace *sevobench*下定义了*single_algorithm_parameter*作为主模版类

```C++
namespace sevobench {


template<std::uint64_t Alg_HashName,std::floating_point T>
struct single_algorithm_parameter;

}
```

+ 模版类的第一个参数*std::uint64_t Alg_HashName*代表着算法函数的*Hash*值,因此每个算法函数都有一个唯一的*Alg_HashName*值,而且不同算法函数之间*Alg_HashName*值不同,从而偏特化*single_algorithm*类,从而实现不同算法由主模版*single_algorithm*类去表达.这样设计的目的是方便将不同的算法在测试集上测试.

+ 其余模版参数的意义与优化算法的模版函数参数意义相同.



&emsp;在*single_algorithm.hpp*头文件中也定义了 
*REGISTER_ALGORITHM(func_name,func,parameter)*
宏表达式


&emsp;*REGISTER_ALGORITHM(func_name,func,parameter)* 让算法函数*func*成为偏特化的*single_algorithm*类,并生成
*func_name*变量,指代算法函数*func*的*Hash*值,偏特化之后的*single_algorithm*类具有一个静态成员字符串数组变量*name*(用于指代算法名字为*func_name*)和一个[运算符重载函数](https://en.cppreference.com/w/cpp/language/operators) *inline auto operator()(F &&f, T l, T r)*,偏特化之后的*single_algorithm*类可以看成一个[函数对象](https://en.cppreference.com/w/cpp/utility/functional
).同理让*parameter*生成一个偏特化之后*single_algorithm_parameter*类,并且继承*parameter*类,详细细节请参考*single_algorithm.hpp*中关于宏表达式的实现代码.


&emsp;下面给出示例,展示用户如何实现一个算法,并进行注册

&emsp;首先包含必要的头文件
```C++
#include"single_algorithm.hpp"


/*
其余实现算法所需要的头文件
*/
```

&emsp;用户实现算法函数*Alogorithm*
```C++


template<std::floating_point T>
struct Algorithm_Paramater {
    /*定义算法所需要的参数*/
};

template<int Dim,int Pop_Size,int Max,bool Memory_Flag,typename F,std::floating_point T,typename Parameter_Type=Algorithm_Paramater<T>>
inline auto Algorithm(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept  {

/*
 声明过程省略
 */

/*
  Max_Itreation取决于模版变量Max,Max是最大迭代次数还是适应度最大计算次数则取决于模版变量B
  */
 for (int iteration=0; iteration < Max_Itreation; iteration++) {
    /*
    中间过程省略
    */
    if constexpr (Memory_Flag) {
            convergence_curve[2*iteration] =Current_FES;//记录当前迭代之后的问题适应度计算次数
            convergence_curve[2*iteration+1] = Current_Best_Fit;//记录当前迭代之后的最佳适应度值
        }
 }


if constexpr(Memory_Flag)
        return std::make_tuple(best_pos, best_fit,convergence_curve);
if constexpr(!Memory_Flag)
    return std::make_pair(best_pos, best_fit);
}
```

&emsp;最后对函数进行注册

```C++
REGISTER_ALGORITHM(Algorithm_Name,Algorithm,Algorithm_Paramater)
```

&emsp;这样用户通过实现*Algorithm*函数,得到偏特化的一个*single_algorithm*类和*single_algorithm_parameter*类,其中偏特化后的*single_algorithm_parameter*类继承*Algorithm_Paramater*,同时得到一个变量*Algorithm_Name*去指代*Algorithm*函数的*Hash*值.

&emsp;我们可以使用偏特化之后的类去优化[Rastrigin问题](https://www.sfu.ca/~ssurjano/rastr.html)
```C++
#define ALGORITHM(ALG,Dim,F,L,R) single_algorithm<ALG, Dim, 30,1000*Dim, false>()(F,L,R)

int main() {
    std::cout<<ALGORITHM(Algorithm_Name, 30,rastrigin<30>, -100.0f, 100.0f).second<<"\n";
}

```



### [测试问题集](#sevobench用户指南)

&emsp;想通过算法连续测试一系列问题时,需要构建关于问题集的类.

&emsp;在SEvoBench中,一般声明一个问题集,使用如下形式去声明

```C++
template< int Prob_Index, int Dim,std::floating_point T >
struct  Problem;
```
&emsp;*Problem*类的设计依然使用[模版偏特化技术](https://en.cppreference.com/w/cpp/language/partial_specialization) 

+ 模版参数 *int Prob_Index*:代表问题索引
+ 模版参数 *int Dim*:代表问题维度
+ *std::floating_point T*:浮点数类型

&emsp;使用第一个模版参数*int Prob_Index*去偏特化*Problem*类.

```C++
template<int Dim,std::floating_point T>
struct Problem<1,Dim,T> {
    static constexpr T L=-100; // 问题下界,使用变量L
    static constexpr T U=100; // 问题上界,使用变量U
   
   /*
   使用运算符重载函数T operator()(T *x)
   */

    T operator()(T *x) noexcept {
       return 10*Dim+std::accumulate(_x, _x+Dim, T(0), [](auto l,auto r){
        return l+r*r-10*std::cos(2*std::numbers::pi_v<T>*r);
    });
    }
};

```

&emsp;我们偏特化一个*Problem*类,它的问题索引为1,问题类型为[Rastrigin问题](https://www.sfu.ca/~ssurjano/rastr.html).

### [问题集主模版类](#sevobench用户指南)
&emsp;为了使得不同的问题集有统一的问题集类去表示,SEvoBench使用[模版偏特化技术](https://en.cppreference.com/w/cpp/language/partial_specialization) 

&emsp;对不同的问题测试集,用一个主模板类去表示. 

&emsp;在*single_problem.hpp*头文件中在namespace *sevobench*下声明了*single_problem*作为主模版类

```C++
namespace sevobench {

template <std::uint64_t Prob_HashName, int Prob_Index, int Dim,
          std::floating_point T>
struct single_problem;

}
```
&emsp;对于不同的问题集,有不同的初始化方法,为此也声明一个主模版类用于初始化问题集.
```C++
namespace sevobench {

template <std::uint64_t Prob_HashName, int Dim, std::floating_point T>
struct init_single_problem;

}
```

+ 和算法基类类似,类的第一个模版参数*std::uint64_t Prob_HashName*为问题集的*Hash*值,用于偏特化.

+ 其余模版参数和测试问题集的模版参数意义一致

&emsp;为此提供了两种宏表达式

`REGISTER_PROBLEM(problem_name,problem,size)`

`INIT_PROBLEM(problem,F)`

分别用于偏特化问题集和初始化方法.


&emsp;下面给出示例,展示用户如何实现一个问题集,并进行偏特化

&emsp;首先包含必要的头文件
```C++
#include"single_problem.hpp"


/*
其余实现算法所需要的头文件
*/
```

&emsp;声明和定义问题,定义初始化方法
```C++
template< int Prob_Index, int Dim,std::floating_point T >
Struct  Problem;



template<int Dim,std::floating_point T>
struct Problem<1,Dim,T> {
    static constexpr T L=-100; // 问题下界,使用变量L
    static constexpr T U=100; // 问题上界,使用变量U
    T operator()(T *x) noexcept {
       return 10*Dim+std::accumulate(_x, _x+Dim, T(0), [](auto l,auto r){
        return l+r*r-10*std::cos(2*std::numbers::pi_v<T>*r);
    });
    }
};

template<std::floating_point T,int Dim>
void init_Problem() {
    /*
    用户定义
    */
}
```

&emsp;偏特化问题和初始化方法 
```C++
REGISTER_PROBLEM(Problem_Name,Problem,1) //最后一个参数代表问题集的问题个数
INIT_PROBLEM(Problem,init_Problem)
```

&emsp;这样用户通过实现*Problem*问题集,并对其进行偏特化为一个*single_problem*类,得到一个变量*Problem_Name*去指代*Problem*问题集的*Hash*值.


&emsp;不过一般来说,用户不需要去添加问题集,SEvoBench平台内置了一系列问题集.




### [算法测试函数](#sevobench用户指南)
&emsp;在*parallel_algorithm_benchmark.hpp*头文件中namespace *sevobench*下,我们定义了
```C++
namespace sevobench {

template<std::uint64_t Alg_HashName,std::uint64_t Prob_HashName,bool Memory_Flag,int Dim,int Pop_Size,int Max,int Runs,bool parallel,std::floating_point T,typename Parameter_Type=typename single_algorithm_parameter<Alg_HashName, T>::base_type>
inline auto evo_bench(const Parameter_Type & par=Parameter_Type()) noexcept;

}
 ```
&emsp;对于用户来说只需要关注这几个模版参数

+ *std::uint64_t Alg_HashName*:算法的*Hash*变量值指定特定算法
+ *std::uint64_t Prob_HashName*:问题集的*Hash*变量值指定特定问题
+ *bool Memory_Flag*:是否记录收敛过程,默认值为false
+ *int Dim*:问题集的维度,默认值30
+ *int Pop_Size*:进化算法的种群数,默认值30
+ *int Max*:问题的适应度最大计算次数或者算法迭代次数,默认为1000*Dim
+ *int Runs*:测试重复次数,默认值是30
+ *bool parallel*:是否并行计算,默认是true
+ *std::floating_point T*:浮点数类型,默认是float
+ *typename Parameter_Type*: 指定封装算法参数的类,默认类型为 *REGISTER_ALGORITHM(Algorithm_Name,Algorithm,Algorithm_Paramater)* 中的 *Algorithm_Paramater*

&emsp;测试函数返回类型为*benchmark_result\<std::floating_point T,typename Algorithm,typename Problem>*,该数据结构存储了算法测试结果.

#### [关于*benchmark_result\<std::floating_point T,typename Algorithm,typename Problem>*](#sevobench用户指南)

+ 定义在*parallel_algorithm_benchmark.hpp*头文件中namespace *sevobench::sevobench_detail*

+ 声明如下
```C++
namespace sevobench::sevobench_detail {

template<std::floating_point T,typename Algorithm,typename Problem>
class benchmark_result;

}
```
+ 具有的成员函数

> *auto get_graph_data()* :返回收敛数据*graph_data*,一维的数组形式存储,数据结构为*std::vector\<T\>*.

> *auto get_table_data()* :返回统计数据*table_data*,一维的数组形式存储,数据结构为*std::vector\<T\>*.

> *constexpr auto problem_dim()* :返回测试集问题的维度

> *constexpr auto problem_size()* :返回测试集问题的数目


> *constexpr auto problem_name()* :返回测试集的名字

> *constexpr auto algorithm_name()* :返回测试算法的名字

> *constexpr auto algorithm_max()* : 返回算法的最大迭代次数或者问题适应度值最大计算次数

> *constexpr auto algorithm_flag()* : 返回算法的*Memory_Flag*标记


 > *void write_table(const char \*dir_name)* :将统计数据写入*dir_name*文件夹下,文件以`{Algorithm}_{Populations}_{Problem}_{Dim}_table.txt`命名
 
 该文件的数据组织形式为
```

...

Mean
Std
Best
Time
25% quantile
median
75% quantile
Worst
...

```
文件每一行一个数据,每8个数据作为一组,每一组索引对应着相应问题集的索引.这8个数据分别代表算法测试第I个问题运行N次之后的平均值,方差,最优值,平均时间,较小四分位数,中位数,较大四分位数,最差值.


> *void write_curve(const char \*dir_name)* :将收敛数据写入*dir_name*文件夹下,生成多个文件,每个文件以`{Algorithm}_{Populations}_{Problem}_F{index}_{Dim}_curve.txt`方式命名

该文件的数据组织形式为
```
...
FES
Best_Fit
···

```
&emsp;文本每一行一个数据,依次为*FES*(问题适应度计算次数),*Best_Fit*(当前最优值),依次循环往复.其中{*FES*,*Best_Fit*}组成一个点,算法内部的每次迭代会生成一个点,这个点一共有*Max*(最大迭代次数)个

+ 值得注意的是成员函数*auto get_graph_data()* 和 *auto get_table_data()* 返回的数据类型排布方式与文本数据排布方式相同


## [概念与约束](#sevobench用户指南)
在*common_concept.hpp*中定义了一些通用的概念,用于约束算法函数声明中的模版参数
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

## [API介绍](#sevobench用户指南)

### [关于问题集的API](#sevobench用户指南)

#### [*ShiftFunc*](#sevobench用户指南)

+ 问题集的哈希变量名:*sevobench::ShiftFunc*
+ 问题集类型:经典的13个可扩展维度的单目标函数,每个函数经过随机偏移处理
+ 问题集相关文献链接:https://ieeexplore.ieee.org/document/771163
+ 所在头文件:*classic_problem.hpp single_problem.hpp*

&emsp;在*single_problem.hpp*中*namespace sevobench*下的相关声明
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<ShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*RotateShiftFunc*](#sevobench用户指南)

+ 问题集的哈希变量名:*sevobench::RotateShiftFunc*
+ 问题集类型:经典的13个可扩展维度的单目标函数,每个函数经过随机旋转偏移处理
+ 问题集相关文献链接:https://ieeexplore.ieee.org/document/771163
+ 所在头文件:*classic_problem.hpp single_problem.hpp*

&emsp;在*single_problem.hpp*中*namespace sevobench*下的相关声明
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<RotateShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*AsyShiftFunc*](#sevobench用户指南)

+ 问题集的哈希变量名:*sevobench::AsyShiftFunc*
+ 问题集类型:经典的13个可扩展维度的单目标函数,每个函数经过随机偏移处理,并使用Non-Linear Transformations进行处理,使其函数具有Ill-conditioning、Irregularity、Symmetry-breaking特征.
+ 问题集相关文献链接:https://ieeexplore.ieee.org/document/771163
+ 函数变换处理相关文献:Dim. Hansen, S. Finck, R. Ros, A. Auger, Real-Parameter Black-Box Optimization Benchmarking 2009: Noiseless Functions Definitions, Tech. rep. RR-6829,
INRIA, 2010.
+ 所在头文件:*classic_problem.hpp single_problem.hpp*

&emsp;在*single_problem.hpp*中*namespace sevobench*下的相关声明
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<AsyShiftFunc,Prob_Index,Dim,T>

}
```
<br>
<br>

#### [*AsyRotateShiftFunc*](#sevobench用户指南)

+ 问题集哈希变量:*sevobench::AsyRotateShiftFunc*
+ 问题集类型:经典的13个可扩展维度的单目标函数,每个函数经过随机旋转偏移处理,并使用Non-Linear Transformations进行处理,使其函数具有Ill-conditioning、irregularity、symmetry-breaking特征.
+ 问题集相关文献链接:https://ieeexplore.ieee.org/document/771163
+ 函数变换处理相关文献:Dim. Hansen, S. Finck, R. Ros, A. Auger, Real-Parameter Black-Box Optimization Benchmarking 2009: Noiseless Functions Definitions, Tech. rep. RR-6829,
INRIA, 2010.
+ 所在头文件:*classic_problem.hpp single_problem.hpp*

&emsp;在*single_problem.hpp*中*namespace sevobench*下的相关声明
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<AsyRotateShiftFunc,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*CEC2020*](#sevobench用户指南)
+ 问题集哈希变量:*sevobench::CEC2020*
+ 问题集类型:IEEE CEC 2020测试集的10个单目标函数.原始测试集对函数维度进行限制,平台中的测试集支持维度可扩展.
+ 问题集相关链接:https://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC2020/CEC2020-2.htm
+ 所在头文件:*cec2020.hpp single_problem.hpp*

&emsp;在*single_problem.hpp*中*namespace sevobench*下的相关声明
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<CEC2020,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### [*CEC2022*](#sevobench用户指南)
+ 问题集哈希变量:*sevobench::CEC2022*
+ 问题集类型:IEEE CEC 2022测试集的12个单目标函数.原始测试集对函数维度进行限制,平台中的测试集支持维度可扩展.
+ 问题集相关链接:https://www3.ntu.edu.sg/home/epnsugan/index_files/CEC2022/CEC2022.htm
+ 所在头文件:*cec2022.hpp single_problem.hpp*


&emsp;在*single_problem.hpp*中*namespace sevobench*下的相关声明
```C++
namespace sevobench {

template<int Prob_Index,int Dim,std::floating_point T>
struct single_problem<CEC2022,Prob_Index,Dim,T>;

}
```
<br>
<br>

#### 相关说明

&emsp;由于每次问题集的初始化过程,都会随机生成一份相关测试数据,所以比较不同算法在同一测试集上的性能,请在问题集的同一初始化下进行比较.

### [关于算法的API](#sevobench用户指南)

#### [差分进化算法](#sevobench用户指南)
+ 算法哈希变量名:*sevobench::DE*
+ 算法相关文献链接:https://link.springer.com/article/10.1023/A:1008202821328
+ 所在头文件:*de.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=de_parameter<T>>
inline auto de_optimize(G &&positions,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=de_parameter<T>>
inline auto de(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```
+ 参数类*de_parameter* :*f*为weight coefficient,*cr*为crossover probability
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

#### [自适应差分进化算法](#sevobench用户指南)
+ 算法哈希变量名:*sevobench::JADE*
+ 算法相关文献链接:https://ieeexplore.ieee.org/document/5208221
+ 所在头文件:*jade.hpp*
+ 函数声明
```C++
namespace sevobench {

    template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=jade_parameter<T>>
inline auto jade_optimize(G &&positions,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type());

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=jade_parameter<T>>
inline auto jade(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;
 
}
```
+ 参数类 *jade_parameter*:参数与论文一致
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

+ 算法哈希变量名:*sevobench::SHADE*
+ 算法相关文献链接:https://ieeexplore.ieee.org/document/6557555  https://ieeexplore.ieee.org/abstract/document/6900380 
+ 所在头文件:*template_shade.hpp*
+ 函数声明
```C++
namespace sevobench {

   template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=shade_parameter<T>>
inline auto shade(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept ;

}
```

<br>
<br>

+ 算法哈希变量名:*sevobench::LSHADE*
+ 算法相关文献链接:https://ieeexplore.ieee.org/abstract/document/6900380 
+ 所在头文件:*template_shade.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size=18*Dim,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=lshade_parameter<T>>
inline auto lshade(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ 算法参数 *shade_parameter* *lshade_parameter* :与论文参数一致

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

#### [粒子群优化算法](#sevobench用户指南)
+ 算法哈希变量名:*sevobench::PSO*
+ 算法相关文献链接:https://ieeexplore.ieee.org/document/699146
+ 所在头文件:*pso.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim, int Pop_Size = 30, int Max = 1000*Dim, bool Memory_Flag = false, typename G,typename F, std::floating_point T,typename Parameter_Type=pso_parameter<T>>
inline auto pso_optimize(G &&positions,F &&f, T left_bound, T right_bound,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim, int Pop_Size = 30, int Max = 1000*Dim, bool Memory_Flag = false,typename F, std::floating_point T,typename Parameter_Type=pso_parameter<T>>
inline auto pso(F &&f, T l, T r,const Parameter_Type& pt=Parameter_Type()) noexcept ;


}
```

+ 算法参数*pso_parameter*
```C++
namespace sevobench {

template<std::floating_point T>
struct pso_parameter {
    T c1 = T(2);
    T c2 = T(2);
    T w_min = T(0.4);
    T w_max = T(0.9);
    T v_max_ratio =  T(0.2);
};

}
```

<br>
<br>

+ 算法哈希变量名:*sevobench::SPSO2007*
+ 算法相关文献链接:https://ieeexplore.ieee.org/document/4223164 https://hal.science/hal-00764996
+ 所在头文件:*template_pso.hpp*
+ 函数声明
```C++
namespace {

template<int Dim,int Pop_Size=30,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=spso2007_parameter<T>>
inline auto spso2007(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

<br>
<br>

+ 算法哈希变量名:*sevobench::SPSO2011*
+ 算法相关文献链接:https://ieeexplore.ieee.org/document/6557848 
https://hal.science/hal-00764996
+ 所在头文件:*template_pso.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size=30,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=spso2011_parameter<T>>
inline auto spso2011(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept

}
```

+ 算法参数*spso2007_parameter*  *spso2011_parameter* :参数与 https://hal.science/hal-00764996 文献参数一致

```C++
namespace sevobench {

template<std::floating_point T>
struct spso2007_parameter {
    T w=T(1)/T(2*std::numbers::ln2_v<T>);
    T c=T(0.5)+std::numbers::ln2_v<T>;
};

template<std::floating_point T>
struct spso2011_parameter {
    T w=T(1)/T(2*std::numbers::ln2_v<T>);
    T c=T(0.5)+std::numbers::ln2_v<T>;
};

}
```
<br>
<br>

#### [用于解决大规模优化的粒子群优化算法变种](#sevobench用户指南)
+ 算法哈希变量名:*sevobench::CSO*
+ 算法相关文献链接:https://ieeexplore.ieee.org/document/6819057
+ 所在头文件:*cso.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size=250,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=cso_parameter<T>>
inline auto cso_optimize(G &&positions,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim,int Pop_Size=250,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=cso_parameter<T>>
inline auto cso(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ 算法参数  *cso_parameter*:参数与论文参数一致
```C++
namespace sevobench {

template<std::floating_point T>
struct cso_parameter {
    T fai=T(0.1);
};

}
```
<br>
<br>



+ 算法哈希变量名:*sevobench::SLPSO*
+ 算法相关文献链接:https://www.sciencedirect.com/science/article/pii/S0020025514008366
+ 所在头文件:*slpso.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=slpso_parameter<T>>
inline auto slpso_optimize(G &&positions,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim,int Pop_Size,int Max,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=slpso_parameter<T>>
inline auto slpso(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept ;

}
```
+ 算法参数  *slpso_parameter*:参数与论文参数一致
```C++
namespace sevobench {

template<std::floating_point T>
struct slpso_parameter {
    T beta=T(0.01);
    T alpha=T(0.5);
}
;

}

```
<br>
<br>



#### [生物启发式算法](#sevobench用户指南)
+ 人工蜂群算法
+ 算法哈希变量名:*sevobench::ABC*
+ 相关文献

>D. Karaboga, AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL OPTIMIZATION,TECHNICAL REPORT-TR06, Erciyes University, Engineering Faculty, Computer Engineering Department 2005.

>D. Karaboga, B. Basturk, A powerful and Efficient Algorithm for Numerical Function Optimization: Artificial Bee Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39, Issue:3,pp:459-171, November 2007,ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x 

>D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony (ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, Pages 687-697. 

>D. Karaboga, B. Akay, A Comparative Study of Artificial Bee Colony Algorithm,  Applied Mathematics and Computation, 214, 108-132, 2009. 

+ 所在头文件:*abc.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=abc_parameter<T>>
inline auto abc_optimize(G&& sol,F &&f,T d,T w,const Parameter_Type &pt=Parameter_Type()) noexcept ;

template<int Dim,int Pop_Size=100,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=abc_parameter<T>>
inline auto abc(F &&f,T d,T w,const Parameter_Type &pt=Parameter_Type()) noexcept;

}
```

+ 算法参数  *abc_parameter*
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


#### [局部搜索优化算法](#sevobench用户指南)


+ 随机搜索算法
+ 算法哈希变量名:*sevobench::RandomSearch*
+ 所在头文件:random_search.hpp
+ 函数声明
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

+ (1+1)进化策略
+ 算法哈希变量名:*sevobench::ES*
+ 算法相关文献链接: https://link.springer.com/article/10.1023/A:1015059928466
+ 所在头文件:*es.hpp*
+ 函数声明
```C++
namespace sevobench {

template<int Dim,int Pop_Size=1,int Max=1000*Dim,bool Memory_Flag=false,typename G,typename F,std::floating_point T,typename Parameter_Type=es_parameter<T>>
inline auto es_optimize(G&& sol,F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

template<int Dim,int Pop_Size=1,int Max=1000*Dim,bool Memory_Flag=false,typename F,std::floating_point T,typename Parameter_Type=es_parameter<T>>
inline auto es(F&& function,T l,T r,const Parameter_Type& pt=Parameter_Type()) noexcept;

}
```

+ 算法参数  *es_parameter*
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



## [示例](#sevobench用户指南)
为了方便用户理解,在*example*文件夹下给出了五个示例
+ *example/quickstart*:展示如何快速开始使用SEvoBench框架
+ *example/coco_test*:展示如何将SEvoBench的算法集成到COCO平台的BBOB测试集上
+ *example/cmaes_test*: 展示如何将GitHub上的CMA-ES算法开源代码集成到SEvoBench框架上
+ *example/parallel_test*: 展示SEvoBench的并行测试性能
+ *example/add_problem*:展示如何添加问题到SEvoBench框架

