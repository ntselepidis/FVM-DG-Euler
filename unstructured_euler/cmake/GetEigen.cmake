include(ExternalProject)

add_library(Eigen INTERFACE)

ExternalProject_Add(
    EigenProject
    URL https://github.com/eigenteam/eigen-git-mirror/archive/3.3.7.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

add_dependencies(Eigen EigenProject)
target_include_directories(Eigen
  SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/Eigen
)
