#ifndef HYPSYS1D_RECONSTRUCTION_HPP
#define HYPSYS1D_RECONSTRUCTION_HPP

#include <Eigen/Dense>
#include <cmath>
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/limiters.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>


class PWConstantReconstruction {
  public:
    void set(const Eigen::MatrixXd &u) const {
        up.resize(u.rows(),u.cols());
        up = u;
    }

    /// Compute the left and right trace at the interface i + 1/2.
    /** Note: This API is agnostic to the number of cell-averages required
     *        by the method. Therefore, reconstructions with different stencil
     *        sizes can implement this API; and this call can be used in parts
     *        of the code that do not need to know about the details of the
     *        reconstruction.
     */
    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    operator()(int i) const {
        return (*this)(up.col(i), up.col(i+1));
    }

    /// Compute the left and right trace at the interface.
    /** Piecewise constant reconstruction of the left and right trace only
     *  requires the cell-average to the left and right of the interface.
     *
     *  Note: Compared to the other overload this reduces the assumption on
     *        how the cell-averages are stored. This is useful when testing and
     *        generally makes the function useful in more situations.
     */
    inline
    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    operator()(Eigen::VectorXd ua, Eigen::VectorXd ub) const {
        return {std::move(ua), std::move(ub)};
    }

private:
  mutable Eigen::MatrixXd up;
};

template <class SlopeLimiter>
class PWLinearReconstruction {
  public:
    explicit PWLinearReconstruction(const SlopeLimiter &slope_limiter)
        : slope_limiter(slope_limiter) {}

    void set(const Eigen::MatrixXd &u) const {
        up.resize(u.rows(),u.cols());
        up = u;
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd> 
    operator()(int i) const {
        return (*this)(up.col(i-1), up.col(i), up.col(i+1), up.col(i+2));
    }

    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    operator()(Eigen::VectorXd ua, 
               Eigen::VectorXd ub, 
               Eigen::VectorXd uc, 
               Eigen::VectorXd ud) const {

        Eigen::VectorXd sL = ub - ua;
        Eigen::VectorXd sM = uc - ub;
        Eigen::VectorXd sR = ud - uc;
        Eigen::VectorXd uL(3);
        Eigen::VectorXd uR(3);

        for (int i=0; i<3; i++) {
            uL(i) = ub(i) + 0.5 * slope_limiter(sL(i), sM(i));
            uR(i) = uc(i) - 0.5 * slope_limiter(sM(i), sR(i));            
        }

        return {std::move(uL), std::move(uR)};
    }

  private:
    SlopeLimiter slope_limiter;
    mutable Eigen::MatrixXd up;
};

#endif // HYPSYS1D_RATE_OF_CHANGE_HPP
