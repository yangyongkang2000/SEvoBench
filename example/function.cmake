macro(common arg1)
set_target_properties(${arg1} PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED ON)
target_include_directories(${arg1} PRIVATE ../../include ../include ../../../include)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


if(CMAKE_CONFIGURATION_TYPES)
    target_compile_options(${arg1} PRIVATE
            $<$<CONFIG:Release>:
            $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /arch:AVX2 /DINSTRSET=8 /F 8388608 /W4>

            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<AND:$<CXX_COMPILER_ID:Clang>,$<NOT:$<CXX_COMPILER_ID:AppleClang>>>>:
            -ffast-math
            -march=native
            -Wall
            -Wextra
            -Wpedantic
            >

            $<$<CXX_COMPILER_ID:AppleClang>:
            -ffast-math
            -mavx2
            -mfma
            -Wall
            -Wextra
            -Wpedantic
            >
            >
    )
else()
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${arg1} PRIVATE
            $<$<CONFIG:Release>:
            $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /arch:AVX2 /DINSTRSET=8 /F 8388608 /W4>

            $<$<OR:$<CXX_COMPILER_ID:GNU>,$<AND:$<CXX_COMPILER_ID:Clang>,$<NOT:$<CXX_COMPILER_ID:AppleClang>>>>:
            -ffast-math
            -march=native
            -Wall
            -Wextra
            -Wpedantic
            >

            $<$<CXX_COMPILER_ID:AppleClang>:
            -ffast-math
            -mavx2
            -mfma
            -Wall
            -Wextra
            -Wpedantic
            >
            >
    )
endif()
endif()


target_link_libraries(${arg1} PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:pthread>
)
endmacro()