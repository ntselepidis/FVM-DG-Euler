#pragma once

#include <fmt/format.h>

#include "euler.hpp"
#include "mesh.hpp"

class CFLCondition {
  public:
    explicit CFLCondition(const Mesh &mesh) : dx(mesh.getMinimumInradius()) {}

    double operator()(const Eigen::MatrixXd &U) const {
        // compute the cfl condition here.
        // you can use `assert_valid_timestep` to check if
        // the computed value is valid.
        return 0.001;
    }

    void assert_valid_timestep(double dt_cfl) const {
        if (dt_cfl <= 0.0 || !std::isfinite(dt_cfl)) {
            throw std::runtime_error(
                fmt::format("Non-positive timestep: dt = {:.3e}", dt_cfl));
        }
    }

  private:
    double dx;
    double cfl_number = 0.45;
};