#include <ancse/cfl_condition.hpp>

#include <Eigen/Dense>

// FVM
FVM_CFLCondition::FVM_CFLCondition(const Grid &grid,
                                   const std::shared_ptr<Model> &model,
                                   double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double FVM_CFLCondition::operator()(const Eigen::MatrixXd &u) const {

  auto n_cells = grid.n_cells;
  auto n_ghost = grid.n_ghost;

  double a_max = 0.0;
  for (int i = n_ghost; i < n_cells - n_ghost; ++i) {
    a_max = std::max(a_max, model->max_eigenvalue(u.col(i)));
  }

  return cfl_number * grid.dx / a_max;
}

// DG Method
DGM_CFLCondition::DGM_CFLCondition(const Grid &grid,
                                   const std::shared_ptr<Model> &model,
                                   const DGHandler &dg_handler,
                                   double cfl_number)
    : grid(grid), model(model), dg_handler(dg_handler), cfl_number(cfl_number) {
}

double DGM_CFLCondition::operator()(const Eigen::MatrixXd &u) const {

  // Extract scaled 0-th order coeffs
  auto u_avg = dg_handler.build_cell_avg(u);

  auto n_cells = grid.n_cells;
  auto n_ghost = grid.n_ghost;

  double a_max = 0.0;
  for (int i = n_ghost; i < n_cells - n_ghost; ++i) {
    a_max = std::max(a_max, model->max_eigenvalue(u_avg.col(i)));
  }

  return cfl_number * grid.dx / a_max;
}

/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const std::shared_ptr<Model> &model,
                   double cfl_number) {
  return std::make_shared<FVM_CFLCondition>(grid, model, cfl_number);
}

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler, double cfl_number) {
  return std::make_shared<DGM_CFLCondition>(grid, model, dg_handler,
                                            cfl_number);
}
