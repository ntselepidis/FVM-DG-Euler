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
    // implement the loop for DG numerical flux term.
}

/// DG volume integral term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_volume_integral(Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
    // implement the loop for DG volume integral.
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
    // Register the other numerical fluxes.

    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}",
                    std::string(config["flux"])));
}

#undef REGISTER_NUMERICAL_FLUX
