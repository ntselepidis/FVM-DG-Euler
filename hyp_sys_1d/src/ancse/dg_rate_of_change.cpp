#include <ancse/dg_rate_of_change.hpp>

#include <Eigen/Dense>
#include <ancse/config.hpp>
#include <ancse/polynomial_basis.hpp>
#include <ancse/dg_handler.hpp>
#include <ancse/numerical_flux.hpp>
#include <fmt/format.h>


/// DG numerical flux term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_numerical_flux (Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
    // Implements the loop for DG numerical flux term.
    // Input:
    //   - u0: ( m x (p+1) ) x N
    // Output:
    //   - dudt ( initialized ): ( m x (p+1) ) x N
    // Note:
    //   - ( m x (p+1) ) x N => ( n_vars x n_coeff ) x n_cells

    // Helpers
    const int n_coeff = 1 + poly_basis.get_degree();
    const int n_cells = grid.n_cells;
    const int n_ghost = grid.n_ghost;
    const int n_vars = model->get_nvars();

    const Eigen::VectorXd phiL = poly_basis( 1.0 );
    const Eigen::VectorXd phiR = poly_basis( 0.0 );
    Eigen::VectorXd uL, uR;
    Eigen::VectorXd fL, fR = Eigen::VectorXd::Zero(n_vars);

    for (int j = n_ghost-1; j < n_cells-n_ghost; j++)
    {
        // Trace values at -+(j+1/2)
        uL = dg_handler.build_sol( u0.col( j ), 1.0 );
        uR = dg_handler.build_sol( u0.col(j+1), 0.0 );

        // Fluxes at (j-1/2) and (j+1/2)
        fL = fR; // Reuse from previous iteration flux at (j-1/2)
        fR = numerical_flux(uL, uR); //   Compute flux at (j+1/2)

        for (int i = 0; i < n_vars; i++)
        {
            dudt.col(j).segment(i*n_coeff, n_coeff)
                    += fL(i)*phiR - fR(i)*phiL;
        }
    }
}

/// DG volume integral term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_volume_integral(Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
    // Implements the loop for DG volume integral.
    // Input:
    //   - u0: ( m x (p+1) ) x N
    // Output:
    //   - dudt ( updated ): ( m x (p+1) ) x N
    // Note:
    //   - ( m x (p+1) ) x N => ( n_vars x n_coeff ) x n_cells

    // Helpers
    const int n_coeff = 1 + poly_basis.get_degree();
    const int n_quad_points = static_cast<int>(quad_points.size());
    const int n_cells = grid.n_cells;
    const int n_ghost = grid.n_ghost;
    const int n_vars = model->get_nvars();

    // Eval basis derivate wrt x (dphi) for all quad points
    Eigen::MatrixXd dphi(n_coeff, n_quad_points);
    for (int k = 0; k < n_quad_points; k++)
    {
        dphi.col(k) = poly_basis.deriv( quad_points(k) );
    }

    // Loop over all cells (except for ghosts)
    for (int j = n_ghost; j < n_cells-n_ghost; j++)
    {
        // Loop over all quadrature points
        for (int k = 0; k < n_quad_points; k++)
        {
            // Build sol u on quad point k ( n_vars-vector )
            Eigen::VectorXd u_k = dg_handler.build_sol(
                         u0.col(j), quad_points(k) );

            // Compute flux on quad point k ( n_vars-vector )
            Eigen::VectorXd f_k = model->flux( u_k );

            // Loop over all components of flux vector
            for (int i = 0; i < n_vars; i++)
            {
                dudt.col(j).segment(i*n_coeff, n_coeff)
                        += quad_weights(k)*f_k(i)*dphi.col(k);
            }
        }
    }
}


#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)          \
    if (config["flux"] == (token)) {                            \
        return std::make_shared< DGRateOfChange<FluxType> >(    \
            grid, model, flux, poly_basis, dg_handler);                     \
    }


std::shared_ptr<RateOfChange> make_dg_rate_of_change(
    const nlohmann::json &config,
    const Grid &grid,
    const std::shared_ptr<Model> &model,
    const PolynomialBasis &poly_basis,
    const DGHandler &dg_handler,
    const std::shared_ptr<SimulationTime> &simulation_time)
{
    // Register the numerical fluxes.
    REGISTER_NUMERICAL_FLUX("central_flux", CentralFlux, CentralFlux(model))
    REGISTER_NUMERICAL_FLUX("rusanov", Rusanov, Rusanov(model))
    REGISTER_NUMERICAL_FLUX("lax_friedrichs",
        LaxFriedrichs, LaxFriedrichs(grid, model, simulation_time))
    REGISTER_NUMERICAL_FLUX("roe", Roe, Roe(model))
    REGISTER_NUMERICAL_FLUX("hll", HLL, HLL(model))
    REGISTER_NUMERICAL_FLUX("hllc", HLLC, HLLC(model))

    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}",
                    std::string(config["flux"])));
}

#undef REGISTER_NUMERICAL_FLUX
