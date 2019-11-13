#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

using EulerState = Eigen::Vector4d;
using EulerStateGradient = Eigen::Matrix<double, 4, 2>;