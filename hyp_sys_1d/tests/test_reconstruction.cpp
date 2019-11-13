#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ancse/reconstruction.hpp>


TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    Eigen::VectorXd ua(3), ub(3);
    ua << 1.0, 0, 1.50;
    ub << 0.1, 0, 0.15;


    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}


