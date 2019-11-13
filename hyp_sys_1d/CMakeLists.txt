# Mandatory information that CMAKE needs
cmake_minimum_required(VERSION 3.12.0)

# Project name
project(HypSys1D)

# Manage dependencies
include("cmake/GetEigen.cmake")
include("cmake/GetFMT.cmake")
include("cmake/GetJSON.cmake")

include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

add_library(hyp_sys_obj OBJECT)
target_compile_features(hyp_sys_obj PUBLIC cxx_std_17)
target_include_directories(hyp_sys_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hyp_sys_obj PUBLIC Eigen)
target_link_libraries(hyp_sys_obj PUBLIC FMT)
target_link_libraries(hyp_sys_obj PUBLIC JSON)
target_link_libraries(hyp_sys_obj PRIVATE warning_flags)
target_compile_definitions(hyp_sys_obj PUBLIC -D_USE_MATH_DEFINES)

if(CLANG_TIDY_EXE)
  set_target_properties(hyp_sys_obj PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()

add_executable(fvm_burgers)
target_link_libraries(fvm_burgers PRIVATE hyp_sys_obj)

add_executable(fvm_shallow_water)
target_link_libraries(fvm_shallow_water PRIVATE hyp_sys_obj)

add_executable(fvm_euler)
target_link_libraries(fvm_euler PRIVATE hyp_sys_obj)

add_executable(dg_burgers)
target_link_libraries(dg_burgers PRIVATE hyp_sys_obj)

add_executable(dg_euler)
target_link_libraries(dg_euler PRIVATE hyp_sys_obj)

# Testing
include(cmake/GetGTest.cmake)

add_executable(unit_tests)
target_link_libraries(unit_tests hyp_sys_obj)
target_link_libraries(unit_tests GTest)

add_subdirectory(src)
add_subdirectory(tests)
