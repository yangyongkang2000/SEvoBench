macro(common arg1)
    set_target_properties(${arg1} PROPERTIES
            CXX_STANDARD 20
            CXX_STANDARD_REQUIRED ON)
    target_include_directories(${arg1} PRIVATE ../../include ../include ../../../include)
    if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
    endif ()


    if (CMAKE_CONFIGURATION_TYPES)
        target_compile_options(${arg1} PRIVATE
                $<$<CONFIG:Release>:
                $<$<CXX_COMPILER_ID:MSVC>:
                /fp:fast
                /arch:AVX2
                /DINSTRSET=8
                /F 8388608
                /W4
                >

                $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
                -ffast-math
                -Wall
                -Wextra
                -Wpedantic
                $<$<AND:$<PLATFORM_ID:Darwin>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>>:
                -mavx2
                -mfma
                >
                $<$<NOT:$<AND:$<PLATFORM_ID:Darwin>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>>>:
                -march=native
                >
                >
                >
        )
    else ()
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            target_compile_options(${arg1} PRIVATE
                    $<$<CONFIG:Release>:
                    $<$<CXX_COMPILER_ID:MSVC>:
                    /fp:fast
                    /arch:AVX2
                    /DINSTRSET=8
                    /F 8388608
                    /W4
                    >

                    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
                    -ffast-math
                    -Wall
                    -Wextra
                    -Wpedantic
                    $<$<AND:$<PLATFORM_ID:Darwin>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>>:
                    -mavx2
                    -mfma
                    >
                    $<$<NOT:$<AND:$<PLATFORM_ID:Darwin>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>>>:
                    -march=native
                    >
                    >
                    >
            )
        endif ()
    endif ()


    target_link_libraries(${arg1} PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:pthread>
    )
endmacro()