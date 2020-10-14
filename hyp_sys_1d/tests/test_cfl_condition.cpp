#include <Eigen/Dense>
#include <ancse/cfl_condition.hpp>
#include <gtest/gtest.h>

// This will check the CFL condition for Burgers equation.
void check_cfl_condition(const CFLCondition &cfl_condition, const Grid &grid,
                         std::shared_ptr<Model> &model, double cfl_number) {
  int n_ghost = grid.n_ghost;
  int n_cells = grid.n_cells;
  int n_vars = model->get_nvars();

  Eigen::MatrixXd u = Eigen::MatrixXd::Constant(n_vars, n_cells, 1.0);
  for (int i = 0; i < n_cells; ++i) {
    u.col(i) *= -i * i;
  }
  double max_abs_u = (n_cells - n_ghost - 1) * (n_cells - n_ghost - 1);

  for (int i = 0; i < n_ghost; ++i) {
    u.col(i) *= 2.0 * max_abs_u;
    u.col(n_cells - n_ghost + i) *= 2.0 * max_abs_u;
  }

  double dt_cfl_approx = cfl_condition(u);
  double dt_cfl_exact = cfl_number * grid.dx / max_abs_u;

  ASSERT_DOUBLE_EQ(dt_cfl_approx, dt_cfl_exact);
}

TEST(FVM_CFLCondition, Example) {
  auto n_ghost = 2;
  auto n_cells = 10 + 2 * n_ghost;
  auto grid = Grid({0.9, 1.0}, n_cells, n_ghost);
  std::shared_ptr<Model> model = std::make_shared<Burgers>();
  double cfl_number = 0.5;

  auto cfl_condition = FVM_CFLCondition(grid, model, cfl_number);
  check_cfl_condition(cfl_condition, grid, model, cfl_number);
}

TEST(DGM_CFLCondition, Example) {
  auto n_ghost = 2;
  auto n_cells = 10 + 2 * n_ghost;
  auto grid = Grid({0.9, 1.0}, n_cells, n_ghost);
  std::shared_ptr<Model> model = std::make_shared<Burgers>();
  double cfl_number = 0.5;
  PolynomialBasis basis(2);
  DGHandler dg_handler(model, basis);

  auto cfl_condition = DGM_CFLCondition(grid, model, dg_handler, cfl_number);
  check_cfl_condition(cfl_condition, grid, model, cfl_number);
}
