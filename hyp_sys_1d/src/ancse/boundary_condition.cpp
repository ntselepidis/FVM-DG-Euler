#include <ancse/boundary_condition.hpp>

#include <fmt/format.h>

std::shared_ptr<BoundaryCondition>
make_boundary_condition(int n_ghost, const std::string &bc_key) {

    if (bc_key == "periodic") {
        return std::make_shared<PeriodicBC>(n_ghost);
    }

    if (bc_key == "outflow") {
        return std::make_shared<OutflowBC>(n_ghost);
    }

    throw std::runtime_error(
        fmt::format("Unknown boundary condition. [{}]", bc_key));
}
