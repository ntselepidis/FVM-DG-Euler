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

TEST(TestPWLinear, Example) {
  auto rc0 = PWLinearReconstruction{MinMod{}};
  auto rc1 = PWLinearReconstruction{SuperBee{}};
  auto rc2 = PWLinearReconstruction{MonotonizedCentral{}};

  Eigen::VectorXd ua(3), ub(3), uc(3), ud(3);
  ua << 1.0, 0, 1.50;
  ub << 0.1, 0, 0.15;
  uc << 0.1, 0, 0.15;
  ud << 1.0, 0, 1.50;

  auto [uL0, uR0] = rc0(ua, ub, uc, ud);
  auto [uL1, uR1] = rc1(ua, ub, uc, ud);
  auto [uL2, uR2] = rc2(ua, ub, uc, ud);
}
