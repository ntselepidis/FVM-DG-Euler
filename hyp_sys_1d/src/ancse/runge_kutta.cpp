#include <ancse/runge_kutta.hpp>

#include <ancse/includes.hpp>
#include <ancse/config.hpp>
#include <fmt/format.h>

#define REGISTER_FVM_RUNGE_KUTTA(token, RKType)                                \
    if (rk_key == token) {                                                     \
        return std::make_shared<RKType>(                                       \
            rate_of_change, boundary_condition, nullptr, n_vars, n_cells);     \
    }

/// make Runge Kutta for FVM
std::shared_ptr<RungeKutta>
make_runge_kutta(const nlohmann::json &config,
                 const std::shared_ptr<RateOfChange> &rate_of_change,
                 const std::shared_ptr<BoundaryCondition> &boundary_condition,
                 int n_vars,
                 int n_cells)
{
    std::string rk_key = config["time_integrator"];

    REGISTER_FVM_RUNGE_KUTTA("forward_euler", ForwardEuler)

    // Register your SSP2 class.

    throw std::runtime_error(
        fmt::format("Unknown time-integrator. [{}]", rk_key));
}

#undef REGISTER_FVM_RUNGE_KUTTA

#define REGISTER_DG_RUNGE_KUTTA(token, RKType)                                 \
    if (rk_key == token) {                                                     \
        return std::make_shared<RKType>(                                       \
            rate_of_change, boundary_condition, dg_limiting, n_vars, n_cells); \
    }

/// make Runge Kutta for DG
std::shared_ptr<RungeKutta>
make_runge_kutta(const nlohmann::json &config,
                 const std::shared_ptr<RateOfChange> &rate_of_change,
                 const std::shared_ptr<BoundaryCondition> &boundary_condition,
                 const std::shared_ptr<Limiting> &dg_limiting,
                 int n_vars,
                 int n_cells)
{
    std::string rk_key = config["time_integrator"];

    REGISTER_DG_RUNGE_KUTTA("forward_euler", ForwardEuler)

    // Register your SSP2 class.

    throw std::runtime_error(
        fmt::format("Unknown time-integrator. [{}]", rk_key));
}

#undef REGISTER_DG_RUNGE_KUTTA
