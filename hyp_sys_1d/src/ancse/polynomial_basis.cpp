#include "../../include/ancse/polynomial_basis.hpp"


/// Computes the Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: operator() (double xi) const
{
	Eigen::VectorXd poly(p+1);

	// normalized shifted Legendre polynomials
	for (int k = 0; k < p+1; k++) {
		poly(k) = std::sqrt(2*k+1) * std::legendre(k,2*xi-1);
	}

    return scaling_factor*poly;
}

/// Computes the derivative of Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: deriv (double xi) const
{
	Eigen::VectorXd dpoly(p+1);

	// derivatives of normalized shifted Legendre polynomials
	if ( p == 0 ) {
		dpoly(0) = 0.0;
	} else if ( p == 1 ) {
		dpoly(0) = 0.0;
		dpoly(1) = std::sqrt(3) * 2.0;
	} else if ( p == 2 ) {
		dpoly(0) = 0.0;
		dpoly(1) = std::sqrt(3) * 2.0;
		dpoly(2) = std::sqrt(5) * (12*xi - 6);
	} else {
		throw std::runtime_error(
			"Usupported Legendre polynomial degree (>2).");
	}

    return scaling_factor*dpoly;
}
