#include <gtest/gtest.h>
#include "types.hpp"
#include "euler.hpp"
#include "hllc.hpp"

TEST(TestHLLCFlux, consistency) {

	const int n_vars = 4;
	const int n_cols = 4;

	Eigen::MatrixXd u(n_vars, n_cols);

    u.col(0) << 1.0,      0.0,     0.0, 1.5;
    u.col(1) << 0.125,    0.25,    0.0, 0.40;
    u.col(2) << 0.11432, -0.11432, 0.0, 0.26432;

    const double rho = 1.0;
    const double  vx = 2.0;
    const double  vy = 3.0;
    const double   p = 4.0;

	u.col(3) << rho, rho*vx, rho*vy, euler::energyFromPVars(rho, vx, vy, p);

	double TOL = 1E-10;
	for (int k = 0; k < u.cols(); k++) {
	    ASSERT_LE( (euler::flux(u.col(k)) - hllc(u.col(k), u.col(k))).norm(), TOL);
	}
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
