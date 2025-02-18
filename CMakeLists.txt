cmake_minimum_required(VERSION 3.15)

project(SEvoBench VERSION 2.0.0 LANGUAGES CXX DESCRIPTION "SEvoBench: A C++ Framework for Evolutionary Single-Objective Optimization Benchmarking")

option(BUILD_TESTS "Whether to build test cases with ctest" OFF)
option(BUILD_EXAMPLE "Whether to build examples" OFF)

add_library(SEvoBench INTERFACE)

include(GNUInstallDirs)

install(TARGETS SEvoBench
        EXPORT SEvoBenchTargets
       )
       
target_include_directories(SEvoBench
       INTERFACE
       "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
       "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(LIB_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/SEvoBench")

install(EXPORT SEvoBenchTargets
        NAMESPACE SEvoBench::
        DESTINATION "${LIB_INSTALL_DIR}/cmake/"
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SEvoBenchConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/SEvoBenchConfig.cmake"
        INSTALL_DESTINATION "${LIB_INSTALL_DIR}/cmake/"
        PATH_VARS LIB_INSTALL_DIR
)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/SEvoBenchConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/SEvoBenchConfigVersion.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/SEvoBenchConfig.cmake" DESTINATION "${LIB_INSTALL_DIR}/cmake/")

if(BUILD_TESTS)
       add_subdirectory(tests)
endif ()

if(BUILD_EXAMPLE)
       add_subdirectory(example)
endif ()

set(CPACK_PACKAGE_VENDOR "Yongkang Yang")
set(CPACK_PACKAGE_CONTACT "yangyongkang@ustl.edu.cn")
include(CPack)
