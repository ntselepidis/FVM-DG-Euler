#include "../../include/ancse/polynomial_basis.hpp"


/// Computes the Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: operator() (double xi) const
{
    return Eigen::VectorXd::Zero(p+1);
}

/// Computes the derivative of Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: deriv (double xi) const
{
    return Eigen::VectorXd::Zero(p+1);
}
