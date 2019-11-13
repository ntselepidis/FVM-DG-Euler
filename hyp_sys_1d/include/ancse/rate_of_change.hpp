#ifndef HYPSYS1D_RATE_OF_CHANGE_HPP
#define HYPSYS1D_RATE_OF_CHANGE_HPP

#include <Eigen/Dense>

/// Interface for rate of change 'du/dt' used in `RungeKutta`.
class RateOfChange {
  public:
    virtual ~RateOfChange() = default;

    /// Set `dudt` to be the rate of change at `u0`.
    /** Given a semi-discrete formulation of the PDE
     *     du/dt = L(u),
     *  store L(u0) into `dudt`.
     */
    virtual void operator()(Eigen::MatrixXd &dudt,
                            const Eigen::MatrixXd &u0) const = 0;
};

#endif // HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
