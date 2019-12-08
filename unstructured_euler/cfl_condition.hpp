#pragma once

#include <fmt/format.h>

#include "euler.hpp"
#include "mesh.hpp"

class CFLCondition {
  public:
    explicit CFLCondition(const Mesh &mesh) : dx(mesh.getMinimumInradius()) {}

    double operator()(const Eigen::MatrixXd &U) const {

        double a_max = 0.0;
        for (int i = 0; i < U.rows(); i++) {
            a_max = std::max( a_max, euler::maxEigenValue( U.row(i) ) ); // TODO: Check implicit conversion from row vector to col vector (i.e. EulerState)
        }

        double dt = cfl_number * dx / a_max;

        assert_valid_timestep( dt );

        return dt;
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
