#include <ancse/cfl_condition.hpp>

#include <Eigen/Dense>

StandardCFLCondition::StandardCFLCondition(const Grid &grid,
                                           const std::shared_ptr<Model> &model,
                                           double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double StandardCFLCondition::operator()(const Eigen::MatrixXd &u) const {

    auto n_cells = grid.n_cells;
    auto n_ghost = grid.n_ghost;

    double a_max = 0.0;
    for (int i = grid.n_ghost; i < n_cells - n_ghost; ++i) {
        a_max = std::max( a_max, model->max_eigenvalue( u.col(i) ) );
    }

    return cfl_number * grid.dx / a_max;
}

/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   double cfl_number) {
    return std::make_shared<StandardCFLCondition>(grid, model, cfl_number);
}

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler,
                   double cfl_number) {
    // implement this 'factory' for your CFL condition.
    return nullptr;
}
