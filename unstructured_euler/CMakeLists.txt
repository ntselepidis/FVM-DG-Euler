cmake_minimum_required(VERSION 3.12.0)

# Project name
project(unstructured_euler)

find_package(Threads REQUIRED)

# Manage dependencies
include("cmake/GetEigen.cmake")
include("cmake/GetIGL.cmake")
include("cmake/GetFMT.cmake")

include("cmake/fast_debug.cmake")
include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

add_library(third_party INTERFACE)
target_link_libraries(third_party INTERFACE ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(third_party INTERFACE IGL Eigen FMT)

# So that Visual Studio will define M_PI
target_compile_definitions(third_party INTERFACE -D_USE_MATH_DEFINES)

# OpenMP
if(HAS_OPENMP)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    find_package(OpenMP COMPONENTS CXX REQUIRED)

    target_compile_definitions(third_party INTERFACE ANCSE_HAS_OPENMP=1)
    target_link_libraries(third_party INTERFACE OpenMP::OpenMP_CXX)
endif()

# Sampling performance profiler
if(HAS_GPERF)
    target_link_libraries(third_party INTERFACE -lprofiler)
endif()

add_library(unstructured_euler INTERFACE)
target_link_libraries(unstructured_euler INTERFACE warning_flags)
target_compile_features(unstructured_euler INTERFACE cxx_std_17)
target_link_libraries(unstructured_euler INTERFACE third_party)

# Where to find the input data files.
target_compile_definitions(unstructured_euler INTERFACE -DANCSE_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/mesh")

# Naca airfoil experiment
add_executable(naca_airfoil naca_airfoil.cpp)
target_link_libraries(naca_airfoil PUBLIC unstructured_euler)

# (Not quite) isentropic vortex experiment
add_executable(vortex vortex.cpp)
target_link_libraries(vortex PUBLIC unstructured_euler)

# Static code analysis
if(CLANG_TIDY_EXE)
    set_target_properties(vortex PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
    set_target_properties(naca_airfoil PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()

# Unit-testing framework
if(HAS_GTEST)
    include(cmake/GetGTest.cmake)

    add_executable(unit_test unit_test.cpp)
    target_link_libraries(unit_test unstructured_euler GTest)
endif()
