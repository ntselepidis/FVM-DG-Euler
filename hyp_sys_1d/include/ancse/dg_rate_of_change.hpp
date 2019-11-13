#ifndef HYPSYS1D_DG_RATE_OF_CHANGE_HPP
#define HYPSYS1D_DG_RATE_OF_CHANGE_HPP

#include <memory>
#include <Eigen/Dense>
#include <iostream>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

#include <ancse/polynomial_basis.hpp>
#include <ancse/dg_handler.hpp>

/// Compute the rate of change due to DG method.
/** The semidiscrete approximation of a PDE using DG is
 *      du_i/dt = - (<F_{i+0.5},\phi^{-}> - <F_{i-0.5},\phi^{+}>)
 *                + \int_K <F, d/dx {\phi}>
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 */
template <class NumericalFlux>
class DGRateOfChange : public RateOfChange {
  public:

    DGRateOfChange (const Grid &grid,
                    const std::shared_ptr<Model> &model,
                    const NumericalFlux &numerical_flux,
                    const PolynomialBasis &poly_basis,
                    const DGHandler &dg_handler)
        : grid (grid),
          model (model),
          numerical_flux (numerical_flux),
          poly_basis (poly_basis),
          dg_handler (dg_handler)
    {
        std::tie(quad_points, quad_weights) = dg_handler.get_quadrature();
    }

    virtual void operator() (Eigen::MatrixXd &dudt,
                             const Eigen::MatrixXd &u0) const override {
        dudt.setZero();
        eval_numerical_flux(dudt, u0);
        eval_volume_integral(dudt, u0);
    }

    void eval_numerical_flux (Eigen::MatrixXd &dudt,
                              const Eigen::MatrixXd &u0) const;

    void eval_volume_integral (Eigen::MatrixXd &dudt,
                               const Eigen::MatrixXd &u0) const;


  private:
    Grid grid;
    std::shared_ptr<Model> model;
    NumericalFlux numerical_flux;
    PolynomialBasis poly_basis;
    DGHandler dg_handler;

    mutable Eigen::VectorXd quad_points, quad_weights;
};

std::shared_ptr<RateOfChange>
make_dg_rate_of_change(const nlohmann::json &config,
                       const Grid &grid,
                       const std::shared_ptr<Model> &model,
                       const PolynomialBasis &poly_basis,
                       const DGHandler &dg_handler,
                       const std::shared_ptr<SimulationTime> &simulation_time);

#endif // HYPSYS1D_DG_RATE_OF_CHANGE_HPP
