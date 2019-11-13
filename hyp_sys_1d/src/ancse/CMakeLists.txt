target_sources(hyp_sys_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/boundary_condition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cfl_condition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/config.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dg_handler.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dg_limiting.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/dg_rate_of_change.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/fvm_rate_of_change.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/grid.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/model.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/polynomial_basis.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/runge_kutta.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/simulation_time.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/snapshot_writer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/time_loop.cpp
)

