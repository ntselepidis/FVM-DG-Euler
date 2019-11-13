#ifndef HYPSYS1D_CFL_CONDITION_HPP
#define HYPSYS1D_CFL_CONDITION_HPP

#include <cmath>
#include <limits>
#include <memory>

#include <Eigen/Dense>

#include <ancse/includes.hpp>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/dg_handler.hpp>

/// Interface for computing CFL restricted timesteps.
class CFLCondition {
  public:
    virtual ~CFLCondition() = default;

    /// Compute the largest time-step satisfying the CFL condition.
    virtual double operator()(const Eigen::MatrixXd &u) const = 0;
};


/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   double cfl_number);

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid,
                   const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler,
                   double cfl_number);

#endif // HYPSYS1D_CFL_CONDITION_HPP
