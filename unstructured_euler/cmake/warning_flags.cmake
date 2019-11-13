add_library(warning_flags INTERFACE)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU|Clang|AppleClang")
  # Clang & GCC
  target_compile_options(warning_flags INTERFACE -Wall)
  target_compile_options(warning_flags INTERFACE -Wextra)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # Intel C++
  # add if using
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # MSVC
  target_compile_options(warning_flags INTERFACE /W3)
endif()
