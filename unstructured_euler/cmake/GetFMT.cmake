add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third_party/fmt-6.0.0)

add_library(FMT INTERFACE)
target_link_libraries(FMT INTERFACE fmt::fmt)
