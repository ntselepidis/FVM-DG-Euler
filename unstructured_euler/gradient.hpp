#pragma once

#include "types.hpp"
#include <Eigen/Dense>

/// Compute the gradients of q_bar.
/** Note: this routine does not care if q are the conserved
 *  or primitive variables.
 *
 * @param [out] dqdx  approximation of dq/dx. Has shape (n_cells, 4).
 * @param [out] dqdy  approximation of dq/dy. Has shape (n_cells, 4).
 * @param       q_bar the cell-averages of `q`. Has shape (n_cells, 4).
 * @param       mesh
 */
void compute_gradients(Eigen::MatrixXd &dqdx,
                       Eigen::MatrixXd &dqdy,
                       const Eigen::MatrixXd &q_bar,
                       const Mesh &mesh) {

    // Compute the gradient of all 4 components of q_bar
}
