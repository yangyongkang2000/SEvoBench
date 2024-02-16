# SEvoBench
## A C++ Framework for Evolutionary Single-Objective Optimization Benchmarking

### Introduction
We introduce SEvoBench, a C++ framework for the benchmarking of evolutionary single-objective optimization. The framework is divided into three main modules, which are used for the construction of optimization algorithms and test suites, and the statistical analysis of test results. Each module is independent and cooperative, and uses template programming techniques to make them generic and easy to extend. For the user’s convenience, the framework also incorporates a series of evolutionary algorithms (EAs) and test suites, which are used for the construction and evaluation of optimization methods.
Moreover, SEvoBench supports parallel testing to accelerate the performance testing of EAs. We aim to achieve the following design goals: 
+	The framework offers efficient implementations of EAs and benchmark suites and enables consistent testing of various algorithms and test suites.
+	By focusing on EAs testing, the framework reduces unnecessary complex design of each module, and maintains simplicity, efficiency, flexibility, reusability, and extensibility of the interfaces.
+	The framework is cross-platform, conforms to the C++20 standard, and has minimal dependency on third-party libraries. By employing template programming techniques extensively, the framework allows the compiler to optimize the program fully at compile time , and only needs the relevant header files to be included, without requiring pre-compiled link libraries, thus simplifying the compilation process.
+	The framework supports parallelization of the testing process, reduces testing time, and enables fast acquisition of test results.

### Future Work
SEvoBench is still under development, and some features and functionalities of it need improvement. Moreover, the number of built-in large-scale optimization algorithms and problems is limited, but we will expand them as the development progresses. 

### Documentation
To learn how to use and for more details about the framework, please refer to the documentation.

+ [Chinese Documentation](./doc/User_Guide_Chinese.md)
+ [English Documentation](./doc/User_Guide.md)

