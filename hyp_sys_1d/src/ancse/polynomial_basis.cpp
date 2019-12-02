#include "../../include/ancse/polynomial_basis.hpp"


/// Computes the Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: operator() (double xi) const
{
	Eigen::VectorXd poly(p+1);

	for (int i = 0; i < p+1; i++)
		poly(i) = std::legendre(i,2*xi-1); // shifted Legendre

    return poly;
}

/// Computes the derivative of Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: deriv (double xi) const
{
	Eigen::VectorXd dpoly(p+1);

	dpoly(0) = 0.0;

	if ( p > 0 ) dpoly(1) = 2.0;
	if ( p > 1 ) dpoly(2) = 12*xi - 6;
	if ( p > 2 ) dpoly(3) = 60*xi*xi - 60*xi + 12;
	if ( p > 3 ) throw std::runtime_error(
		"Usupported Legendre polynomial order.");

    return dpoly;
}
