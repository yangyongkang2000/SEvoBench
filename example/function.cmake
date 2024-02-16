macro(common arg1)
set_target_properties(${arg1} PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED ON)
target_include_directories(${arg1} PRIVATE ../../include ../include)
if(NOT CMAKE_BUILD_TYPE)
    set_target_properties(${arg1} PROPERTIES CMAKE_BUILD_TYPE Release)
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
set_target_properties(${arg1} PROPERTIES CMAKE_CONFIGURATION_TYPES Release)
endif ()
target_compile_options(${arg1} PRIVATE $<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:MSVC>:/fp:fast /arch:AVX2 /F 8388608>> $<$<CONFIG:Release>:
$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-ffast-math -march=native -Wall -Wextra -Wpedantic>>)
target_link_libraries(${arg1} PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:pthread>
)
endmacro()