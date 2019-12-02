#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <iostream>
#include <ancse/polynomial_basis.hpp>

TEST(TestLegendre, Example) {

	PolynomialBasis B0(0);
	PolynomialBasis B1(1);
	PolynomialBasis B2(2);
	PolynomialBasis B3(3);

	Eigen::VectorXd L0(1);
	Eigen::VectorXd L1(2);
	Eigen::VectorXd L2(3);
	Eigen::VectorXd L3(4);

	Eigen::VectorXd xis = Eigen::VectorXd::LinSpaced(10,0.0,1.0);

	for (int j=0; j<xis.cols(); j++) {
		const double xi = xis(j);
		L0 << 1.0;
		L1 << 1.0, 2*xi-1;
		L2 << 1.0, 2*xi-1, 6*xi*xi-6*xi+1;
		L3 << 1.0, 2*xi-1, 6*xi*xi-6*xi+1, 
			20*xi*xi*xi - 30*xi*xi + 12*xi - 1;
		ASSERT_EQ( B0(xi), L0 );
		ASSERT_EQ( B1(xi), L1 );
		ASSERT_EQ( B2(xi), L2 );
		ASSERT_EQ( B3(xi), L3 );
	}
}

TEST(TestLegendreDeriv, Example) {

	PolynomialBasis B0(0);
	PolynomialBasis B1(1);
	PolynomialBasis B2(2);
	PolynomialBasis B3(3);

	Eigen::VectorXd L0(1);
	Eigen::VectorXd L1(2);
	Eigen::VectorXd L2(3);
	Eigen::VectorXd L3(4);

	Eigen::VectorXd xis = Eigen::VectorXd::LinSpaced(10,0.0,1.0);

	for (int j=0; j<xis.cols(); j++) {
		const double xi = xis(j);
		L0 << 0.0;
		L1 << 0.0, 2.0;
		L2 << 0.0, 2.0, 12*xi-6;
		L3 << 0.0, 2.0, 12*xi-6, 60*xi*xi - 60*xi + 12;
		ASSERT_EQ( B0.deriv(xi), L0 );
		ASSERT_EQ( B1.deriv(xi), L1 );
		ASSERT_EQ( B2.deriv(xi), L2 );
		ASSERT_EQ( B3.deriv(xi), L3 );
	}
}
