target_sources(unit_tests
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_model.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_boundary_condition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_cfl_condition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_grid.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_numerical_flux.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/test_reconstruction.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/unit_tests.cpp
)

