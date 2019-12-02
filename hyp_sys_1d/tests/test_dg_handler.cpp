#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <iostream>
#include <ancse/dg_handler.hpp>

TEST(TestDGHandler, Example) {

	std::shared_ptr<Model> model_euler = std::make_shared<Euler>();

	const int n_vars = model_euler->get_nvars();
	const int n_cells = 10;
	const int deg = 2;
	const int n_coeff = deg + 1;
	
	PolynomialBasis basis(deg);

	DGHandler dg_handler( model_euler, basis );

	Eigen::MatrixXd u(n_vars*n_coeff, n_cells);

	for (int j = 0; j < u.cols(); j++) {
		for (int i = 0; i < u.rows(); i++) {
			u(i,j) = j*u.rows() + i;
		}
	}

	Eigen::MatrixXd u_avg_target(n_vars, n_cells);
	
	for (int j = 0; j < n_cells; j++) {
		for (int i = 0; i < n_vars; i++) {
			u_avg_target(i,j) = (j*u.rows() + i) + 2*i;
		}
	}

	auto u_avg = dg_handler.build_cell_avg(u);

	ASSERT_EQ( u_avg, u_avg_target );
}
