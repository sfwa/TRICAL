cmake_minimum_required(VERSION 3.10)
project(TRICAL VERSION 2.0.0 LANGUAGES C)

# Floating-point type used for the running sums and the fit; see TRICAL.h.
set(TRICAL_WIDE_TYPE "double" CACHE STRING
    "Floating-point type for TRICAL's running sums and fit (double or float)")
option(TRICAL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(TRICAL_SOURCES src/TRICAL.c src/fit.c)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    set(TRICAL_WARNINGS -Wall -Wextra -Wpedantic -Wshadow -Wconversion
        -Wdouble-promotion -Wstrict-prototypes -Wmissing-prototypes
        -Wcast-qual)
    if(TRICAL_WARNINGS_AS_ERRORS)
        list(APPEND TRICAL_WARNINGS -Werror)
    endif()
endif()

function(trical_configure target)
    target_include_directories(${target}
        PUBLIC include
        PRIVATE src)
    set_target_properties(${target} PROPERTIES
        C_STANDARD 99
        C_STANDARD_REQUIRED ON
        C_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON)
    target_compile_options(${target} PRIVATE ${TRICAL_WARNINGS})
    # The library and its users must agree on the layout of the instance.
    if(NOT TRICAL_WIDE_TYPE STREQUAL "double")
        target_compile_definitions(${target}
            PUBLIC TRICAL_WIDE_TYPE=${TRICAL_WIDE_TYPE})
    endif()
    if(UNIX)
        target_link_libraries(${target} PUBLIC m)
    endif()
endfunction()

add_library(TRICAL SHARED ${TRICAL_SOURCES})
trical_configure(TRICAL)
# The public functions carry no export annotations.
set_target_properties(TRICAL PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(TRICALstatic STATIC ${TRICAL_SOURCES})
trical_configure(TRICALstatic)

enable_testing()
add_subdirectory(test EXCLUDE_FROM_ALL)
