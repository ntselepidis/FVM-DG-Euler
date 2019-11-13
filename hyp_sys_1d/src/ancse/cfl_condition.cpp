#include <ancse/cfl_condition.hpp>

#include <Eigen/Dense>


/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   double cfl_number) {
    // implement this 'factory' for your CFL condition.
    return nullptr;
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
