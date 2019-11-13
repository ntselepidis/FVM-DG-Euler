#include <ancse/fvm_rate_of_change.hpp>

#include <Eigen/Dense>
#include <ancse/numerical_flux.hpp>
#include <ancse/reconstruction.hpp>
#include <fmt/format.h>

#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)                         \
    if (config["flux"] == (token)) {                                           \
        return std::make_shared<FVMRateOfChange<FluxType, Reconstruction>>(    \
            grid, model, flux, reconstruction);                                       \
    }

template <class Reconstruction>
std::shared_ptr<RateOfChange>
deduce_numerical_flux(const nlohmann::json &config,
                      const Grid &grid,
                      const std::shared_ptr<Model> &model,
                      const std::shared_ptr<SimulationTime> &simulation_time,
                      const Reconstruction &reconstruction)
{
    REGISTER_NUMERICAL_FLUX("central_flux", CentralFlux, CentralFlux(model))

    // Register the other numerical fluxes.

    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}", std::string(config["flux"])));
}
#undef REGISTER_NUMERICAL_FLUX

#define REGISTER_RECONSTRUCTION(token, reconstruction)                         \
    if (config["reconstruction"] == token) {                                   \
        return deduce_numerical_flux(                                          \
            config, grid, model, simulation_time, reconstruction);                     \
    }

std::shared_ptr<RateOfChange> make_fvm_rate_of_change(
    const nlohmann::json &config,
    const Grid &grid,
    const std::shared_ptr<Model> &model,
    const std::shared_ptr<SimulationTime> &simulation_time)
{
    REGISTER_RECONSTRUCTION("o1", PWConstantReconstruction{})

    // Register piecewise linear reconstructions.

    throw std::runtime_error(fmt::format(
        "Unknown reconstruction. [{}]", std::string(config["reconstruction"])));
}

#undef REGISTER_RECONSTRUCTION
