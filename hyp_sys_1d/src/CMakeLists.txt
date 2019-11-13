add_subdirectory(ancse)

target_sources(fvm_burgers
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fvm_burgers.cpp
)

target_sources(fvm_shallow_water
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fvm_shallow_water.cpp
)

target_sources(fvm_euler
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fvm_euler.cpp
)

target_sources(dg_burgers
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/dg_burgers.cpp
)

target_sources(dg_euler
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/dg_euler.cpp
)
