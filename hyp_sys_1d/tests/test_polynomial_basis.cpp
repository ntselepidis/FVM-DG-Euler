#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ancse/polynomial_basis.hpp>
#include <iostream>

TEST(TestLegendre, Example) {

  const double scal = 2.0;
  PolynomialBasis B0(0, scal), B1(1, scal), B2(2, scal);
  Eigen::VectorXd L0(1), L1(2), L2(3);
  Eigen::VectorXd xis = Eigen::VectorXd::LinSpaced(10, 0.0, 1.0);

  for (int j = 0; j < xis.cols(); j++) {
    const double xi = xis(j);
    L0 << 1.0;
    L1 << 1.0, sqrt(3) * (2 * xi - 1);
    L2 << 1.0, sqrt(3) * (2 * xi - 1), sqrt(5) * (6 * xi * xi - 6 * xi + 1);
    ASSERT_EQ(B0(xi), scal * L0);
    ASSERT_EQ(B1(xi), scal * L1);
    ASSERT_EQ(B2(xi), scal * L2);
  }
}

TEST(TestLegendreDeriv, Example) {

  const double scal = 2.0;
  PolynomialBasis B0(0, scal), B1(1, scal), B2(2, scal);
  Eigen::VectorXd L0(1), L1(2), L2(3);
  Eigen::VectorXd xis = Eigen::VectorXd::LinSpaced(10, 0.0, 1.0);

  for (int j = 0; j < xis.cols(); j++) {
    const double xi = xis(j);
    L0 << 0.0;
    L1 << 0.0, sqrt(3) * 2.0;
    L2 << 0.0, sqrt(3) * 2.0, sqrt(5) * (12 * xi - 6);
    ASSERT_EQ(B0.deriv(xi), scal * L0);
    ASSERT_EQ(B1.deriv(xi), scal * L1);
    ASSERT_EQ(B2.deriv(xi), scal * L2);
  }
}
