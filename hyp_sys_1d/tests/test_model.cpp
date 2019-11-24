#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <iostream>
#include <ancse/model.hpp>

TEST(TestEulerPrimToCons, Example) {

	Euler model_euler;

	const double gamma = model_euler.get_gamma();
	
	const double rho = 1.0;
	const double   v = -2.0;
	const double   p = 1.0;
	const double   E = p/(gamma-1) + 0.5*rho*v*v;

	Eigen::VectorXd u_prim(3); 
	u_prim << rho, v, p;

	Eigen::VectorXd u_cons_target(3); 
	u_cons_target << rho, rho*v, E;

	Eigen::VectorXd u_cons = model_euler.prim_to_cons( u_prim );

	for (int i=0; i<3; i++) {
		ASSERT_DOUBLE_EQ( u_cons(i), u_cons_target(i) );
	}
}

TEST(TestEulerConsToPrim, Example) {

	Euler model_euler;

	const double gamma = model_euler.get_gamma();
	
	const double rho = 1.0;
	const double   v = -2.0;
	const double   p = 1.0;
	const double   E = p/(gamma-1) + 0.5*rho*v*v;

	Eigen::VectorXd u_cons(3);
	u_cons << rho, rho*v, E;

	Eigen::VectorXd u_prim_target(3);
	u_prim_target << rho, v, p;

	Eigen::VectorXd u_prim = model_euler.cons_to_prim( u_cons );

	for (int i=0; i<3; i++) {
		ASSERT_DOUBLE_EQ( u_prim(i), u_prim_target(i) );
	}
}

TEST(TestEulerEigenvalues, Example) {

	Euler model_euler;

	const double gamma = model_euler.get_gamma();

	const double rho = 1.0;
	const double   v = -2.0;
	const double   p = 1.0;
	const double   c = std::sqrt( gamma*p ) / rho;

	Eigen::VectorXd u_prim(3);
	u_prim << rho, v, p;

	Eigen::VectorXd u_cons = model_euler.prim_to_cons( u_prim );

	Eigen::VectorXd eigvals_target(3);
	eigvals_target << v-c, v, v+c;

	Eigen::VectorXd eigvals = model_euler.eigenvalues( u_cons );

	for (int i=0; i<3; i++) {
		ASSERT_DOUBLE_EQ( eigvals(i), eigvals_target(i) );
	}
}

TEST(TestEulerEigenvectors, Example) {
	
	Euler model_euler;

	const double gamma = model_euler.get_gamma();

	const double rho = 1.0;
	const double   v = -2.0;
	const double   p = 1.0;
	const double   c = std::sqrt( gamma*p ) / rho;

	Eigen::VectorXd u_prim(3);
	u_prim << rho, v, p;

	Eigen::VectorXd u_cons = model_euler.prim_to_cons( u_prim );

	const double E = u_cons(2);
	const double H = (E+p)/rho;

	Eigen::MatrixXd eigvecs_target(3,3);
    eigvecs_target.col(0) << 1.0, v-c, H-v*c;
    eigvecs_target.col(1) << 1.0, v, 0.5*v*v;
    eigvecs_target.col(2) << 1.0, v+c, H+v*c;

    Eigen::MatrixXd eigvecs = model_euler.eigenvectors( u_cons );

    for (int j=0; j<3; j++)
    	for (int i=0; i<3; i++)
    		ASSERT_DOUBLE_EQ( eigvecs(i,j), eigvecs_target(i,j) );
}

TEST(TestEulerMaxEigenvalue, Example) {

	Euler model_euler;

	const double gamma = model_euler.get_gamma();

	const double rho = 1.0;
	const double   v = -2.0;
	const double   p = 1.0;
	const double   c = std::sqrt( gamma*p ) / rho;

	Eigen::VectorXd u_prim(3);
	u_prim << rho, v, p;

	Eigen::VectorXd u_cons = model_euler.prim_to_cons( u_prim );

	double max_eigval = model_euler.max_eigenvalue( u_cons );

	ASSERT_DOUBLE_EQ( max_eigval, 
		std::max( std::abs(v+c), std::abs(v-c) ) );
}

TEST(TestEulerFlux, Example) {

	Euler model_euler;

	const double rho = 1.0;
	const double   v = -2.0;
	const double   p = 1.0;

	Eigen::VectorXd u_prim(3);
	u_prim << rho, v, p;

	Eigen::VectorXd u_cons = model_euler.prim_to_cons( u_prim );

	Eigen::VectorXd fl = model_euler.flux( u_cons );

	const double E = u_cons(2);

	Eigen::VectorXd fl_target(3);
	fl_target << rho*v, rho*v*v+p, (E+p)*v;

	for (int i=0; i<3; i++)
		ASSERT_DOUBLE_EQ( fl(i), fl_target(i) );
}
