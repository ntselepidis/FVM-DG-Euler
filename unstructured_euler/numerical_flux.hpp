#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include "gradient.hpp"
#include "hllc.hpp"
#include "mesh.hpp"
#include "slope_limiter.hpp"

// Note: this class will compute the rate of change due to the fluxes.
// Note: the reason we made this a class is that it allows you to allocate
//       buffers, once at the beginning of the simulation. Add these buffers
//       as needed.
class FluxRateOfChange {
  public:
    explicit FluxRateOfChange(int n_cells) {
        // Allocate buffers as needed.
    }

    void operator()(Eigen::MatrixXd &dudt,
                    const Eigen::MatrixXd &u,
                    const Mesh &mesh) const {

        // Compute the rate of change of u.
        // Note: Please use the method `computeFlux` to abstract
        // away the details of computing the flux through a
        // given interface.
        // Note: You can use `assert_valid_flux` to check
        // if what `computeFlux` returns makes any sense.
        // Note: Do not assume `dudt` is filled with zeros.
    }

    void assert_valid_flux(const Mesh &mesh,
                           int i,
                           int k,
                           const EulerState &nF) const {
        // This is mostly for debugging (but also important to check in
        // real simulations!): Make sure our flux contribution is not
        // nan (ie. it is not not a number, ie it is a number)
        if (!euler::isValidFlux(nF)) {
            // clang-format off
            throw std::runtime_error(
                "invalid value detected in numerical flux, " + euler::to_string(nF)
                + "\nat triangle: " + std::to_string(i)
                + "\nedge:        " + std::to_string(k)
                + "\nis_boundary: " + std::to_string(!mesh.isValidNeighbour(i, k)));
            // clang-format on
        }
    }

    /// Compute the flux through the k-th interface of cell i.
    EulerState computeFlux(const Eigen::MatrixXd &U,
                           int i,
                           int k,
                           const Mesh &mesh) const {
        auto boundary_type = mesh.getBoundaryType(i, k);

        if (boundary_type == Mesh::BoundaryType::INTERIOR_EDGE) {
            return computeInteriorFlux(U, i, k, mesh);
        } else {
            if (boundary_type == Mesh::BoundaryType::OUTFLOW_EDGE) {
                return computeOutflowFlux(U, i, k, mesh);
            } else /* boundary_type == Mesh::BoundaryType::WING_EDGE */
            {
                return computeReflectiveFlux(U, i, k, mesh);
            }
        }
    }

    /// Compute the outflow flux through the k-th interface of cell i.
    /** Note: you know that edge k is an outflow edge.
     */
    EulerState computeOutflowFlux(const Eigen::MatrixXd &U,
                                  int i,
                                  int k,
                                  const Mesh &mesh) const {
        // Implement the outflow flux boundary condition.
        return EulerState{};
    }

    /// Compute the reflective boundary flux through the k-th edge of cell i.
    /** Note: you know that edge k is a reflective/wall boundary edge.
     */
    EulerState computeReflectiveFlux(const Eigen::MatrixXd &U,
                                     int i,
                                     int k,
                                     const Mesh &mesh) const {

        // Implement the reflective flux boundary condition.
        return EulerState{};
    }

    /// Compute the flux through the k-th interface of cell i.
    /** Note: This edge is an interior edge, therefore approximate the flux
     * through this edge with the appropriate FVM formulas.
     */
    EulerState computeInteriorFlux(const Eigen::MatrixXd &U,
                                   int i,
                                   int k,
                                   const Mesh &mesh) const {
        // Reconstruct the trace values of U and compute
        // the numerical flux through the k-th interface of
        // cell i.
        return EulerState{};
    }


  private:
    // add any member variables you might need here.
};
