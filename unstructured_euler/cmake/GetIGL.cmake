include(ExternalProject)

ExternalProject_Add(
  IGLProject
  URL https://github.com/libigl/libigl/archive/master.zip
  SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IGL
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

add_library(IGL INTERFACE)
add_dependencies(IGL IGLProject)
target_include_directories(IGL SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/IGL/include)
