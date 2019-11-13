include(ExternalProject)

add_library(JSON INTERFACE)

ExternalProject_Add(
    JSONProject
    URL https://github.com/nlohmann/json/archive/v3.7.0.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/JSON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

add_dependencies(JSON JSONProject)
target_include_directories(JSON
    SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/JSON/include
)
