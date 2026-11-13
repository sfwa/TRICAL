cmake_minimum_required(VERSION 3.10)

# The tests are C++ (googletest); the library itself is C.
enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use an installed googletest if there is one, otherwise download a release.
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    if(CMAKE_VERSION VERSION_LESS 3.14)
        message(FATAL_ERROR
            "googletest was not found. Install it (e.g. libgtest-dev), or "
            "use CMake 3.14 or later, which can download it.")
    endif()
    include(FetchContent)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
        GIT_SHALLOW TRUE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    add_library(GTest::gtest ALIAS gtest)
    add_library(GTest::gtest_main ALIAS gtest_main)
endif()

# Unit tests reach into the library's internals, so they compile the sources
# rather than linking the library.
add_executable(unittest
    ../src/TRICAL.c
    ../src/fit.c
    test_TRICAL.cpp
    test_linalg.cpp
    test_fit.cpp
    test_accuracy.cpp
    test_robustness.cpp)
target_include_directories(unittest PRIVATE ../include ../src)
if(NOT TRICAL_WIDE_TYPE STREQUAL "double")
    target_compile_definitions(unittest PRIVATE
        TRICAL_WIDE_TYPE=${TRICAL_WIDE_TYPE})
endif()
target_link_libraries(unittest GTest::gtest GTest::gtest_main)
if(UNIX)
    target_link_libraries(unittest m)
endif()

add_test(NAME unittest COMMAND unittest)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS unittest)

# Timing of the public API; not run as part of the tests.
add_executable(benchmark benchmark.c)
set_target_properties(benchmark PROPERTIES C_STANDARD 99)
target_link_libraries(benchmark TRICALstatic)
