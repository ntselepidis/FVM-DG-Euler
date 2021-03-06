#ifndef HYPSYS1D_CFL_CONDITION_HPP
#define HYPSYS1D_CFL_CONDITION_HPP

#include <cmath>
#include <limits>
#include <memory>

#include <Eigen/Dense>

#include <ancse/dg_handler.hpp>
#include <ancse/grid.hpp>
#include <ancse/includes.hpp>
#include <ancse/model.hpp>

/// Interface for computing CFL restricted timesteps.
class CFLCondition {
public:
  virtual ~CFLCondition() = default;

  /// Compute the largest time-step satisfying the CFL condition.
  virtual double operator()(const Eigen::MatrixXd &u) const = 0;
};

/// Compute the CFL condition for generic models (FVM).
/** Compute the maximum eigenvalue given all the data, and compute the
 *  CFL restricted time-step accordingly.
 *
 *  Note: this requires looping over the entire grid, and computing a reduction,
 *        i.e. max.
 */
class FVM_CFLCondition : public CFLCondition {
public:
  FVM_CFLCondition(const Grid &grid, const std::shared_ptr<Model> &model,
                   double cfl_number);

  virtual double operator()(const Eigen::MatrixXd &u) const override;

private:
  Grid grid;
  std::shared_ptr<Model> model;
  double cfl_number;
};

/// Compute the CFL condition for generic models (DG Method).
/** Compute the maximum eigenvalue given all the data, and compute the
 *  CFL restricted time-step accordingly.
 *
 *  Note: this requires looping over the entire grid, and computing a reduction,
 *        i.e. max.
 */
class DGM_CFLCondition : public CFLCondition {
public:
  DGM_CFLCondition(const Grid &grid, const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler, double cfl_number);

  virtual double operator()(const Eigen::MatrixXd &u) const override;

private:
  Grid grid;
  std::shared_ptr<Model> model;
  DGHandler dg_handler;
  double cfl_number;
};

/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const std::shared_ptr<Model> &model,
                   double cfl_number);

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const std::shared_ptr<Model> &model,
                   const DGHandler &dg_handler, double cfl_number);

#endif // HYPSYS1D_CFL_CONDITION_HPP
