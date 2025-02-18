cmake_minimum_required(VERSION 3.15)
project(quickstart LANGUAGES CXX)
#set(CMAKE_CXX_STANDARD  20)
include("../function.cmake")
add_executable(quickstart
        quickstart.cpp)
common(quickstart)

