#ifndef HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
#define HYPSYS1D_FVM_RATE_OF_CHANGE_HPP

#include <memory>
#include <Eigen/Dense>

#include <ancse/config.hpp>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

/// Compute the rate of change due to FVM.
/** The semidiscrete approximation of a PDE using FVM is
 *      du_i/dt = - (F_{i+0.5} - F_{i-0.5}) / dx.
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 * @tparam Reconstruction see e.g. `PWConstantReconstruction`.
 */
template <class NumericalFlux, class Reconstruction>
class FVMRateOfChange : public RateOfChange {
  public:
    FVMRateOfChange(const Grid &grid,
                    const std::shared_ptr<Model> &model,
                    const NumericalFlux &numerical_flux,
                    const Reconstruction &reconstruction,
                    const bool reconstruct_in_prim)
        : grid(grid),
          model(model),
          numerical_flux(numerical_flux),
          reconstruction(reconstruction),
          reconstruct_in_prim(reconstruct_in_prim) {}

    virtual void operator()(Eigen::MatrixXd &dudt,
                            const Eigen::MatrixXd &u0) const override {
      auto n_cells = grid.n_cells;
      auto n_ghost = grid.n_ghost;
      auto dx = grid.dx;
      auto n_vars = model->get_nvars();
      Eigen::VectorXd fL, fR = Eigen::VectorXd::Zero(n_vars);
      Eigen::VectorXd uL, uR;

      if (reconstruct_in_prim) {
        Eigen::MatrixXd u0_prim(u0.rows(),u0.cols());
        for (int i = 0; i < u0.cols(); i++) {
          u0_prim.col(i) = model->cons_to_prim( u0.col(i) );
        }
        reconstruction.set( u0_prim );
      } else {
        reconstruction.set( u0 );
      }

      for (int i = n_ghost-1; i < n_cells-n_ghost; i++) {
        std::tie(uL, uR) = reconstruction(i);

        if (reconstruct_in_prim) {
          uL = model->prim_to_cons(uL);
          uR = model->prim_to_cons(uR);
        }

        fL = fR;
        fR = numerical_flux(uL, uR);

        dudt.col(i) = (fL - fR) / dx;
      }
    }

  private:
    Grid grid;
    std::shared_ptr<Model> model;
    NumericalFlux numerical_flux;
    Reconstruction reconstruction;
    bool reconstruct_in_prim;
};

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const nlohmann::json &config,
                        const Grid &grid,
                        const std::shared_ptr<Model> &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
