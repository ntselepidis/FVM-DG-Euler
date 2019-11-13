#pragma once
#include "types.hpp"

// To make things simpler, we hard code GAMMA.
constexpr const double GAMMA = 1.4;

namespace euler {

//! Computes the pressure
double pressure(const EulerState &u) {
    double rho = u[0];
    double m1 = u[1];
    double m2 = u[2];
    double E = u[3];

    double u1 = m1 / rho;
    double u2 = m2 / rho;

    double internalEnergy = E - 0.5 * rho * (u1 * u1 + u2 * u2);
    double p = (GAMMA - 1) * internalEnergy;

    return p;
}

double speedOfSound(const EulerState &u) {
    double rho = u[0];
    double p = pressure(u);
    return std::sqrt(GAMMA * p / rho);
}

//! Computes the flux though an interface perpendicular to the x-axis.
EulerState flux(const EulerState &u) {
    double rho = u[0];
    double m1 = u[1];
    double m2 = u[2];
    double E = u[3];

    double u1 = m1 / rho;

    double p = pressure(u);
    EulerState fluxValue;

    fluxValue[0] = u1 * rho;
    fluxValue[1] = u1 * m1 + p;
    fluxValue[2] = u1 * m2;
    fluxValue[3] = u1 * (E + p);

    return fluxValue;
}

//! Computes the energy
double energyFromPVars(double rho, double vx, double vy, double p) {
    return p / (GAMMA - 1) + 0.5 * rho * (vx * vx + vy * vy);
}

//! Computes the maximum eigen value
double maxEigenValue(const EulerState &u) {

    double rho = u[0];
    double v = std::sqrt(u[1] * u[1] + u[2] * u[2]) / rho;
    return v + speedOfSound(u);
}

double maxEigenValue(const EulerState &leftState,
                     const EulerState &rightState) {

    auto leftEigenValue = maxEigenValue(leftState);
    auto rightEigenValue = maxEigenValue(rightState);

    return std::max(leftEigenValue, rightEigenValue);
}

//! Check that all components of u are real-valued, e.g. not inf or nan.
bool isReal(const EulerState &u) {
    return std::isfinite(u[0]) || std::isfinite(u[1]) || std::isfinite(u[2])
           || std::isfinite(u[3]);
}

//! Does `u` have plausible values?
bool isValidState(const EulerState &u) {
    double rho = u[0];
    double p = pressure(u);
    return isReal(u) && rho > 0.0 && p > 0.0;
}

//! Does the flux `f` have plausible values?
bool isValidFlux(const EulerState &f) { return isReal(f); }

std::string to_string(const EulerState &u) {
    std::stringstream ss;
    ss << "\n{"
       << "\n\trho = " << u[0] << "\n\tmx  = " << u[1] << "\n\tmy  = " << u[2]
       << "\n\tE   = " << u[3] << "\n}\n";

    return ss.str();
}

/// Rotate u into a coordinate system aligned an edge with `normal`.
/** Let n be the normal and tau the tangential. Then this computes
 *
 *  u_rot = (rho, u[1:3].dot(n), u[1:3].dot(tau), E)
 *
 */
EulerState localCoordinates(const EulerState &u,
                            const Eigen::Vector2d &normal) {
    EulerState u_rot;

    u_rot[0] = u[0];

    u_rot[1] = u[1] * normal[0] + u[2] * normal[1];
    u_rot[2] = -u[1] * normal[1] + u[2] * normal[0];

    u_rot[3] = u[3];

    return u_rot;
}

/// Rotate u back from the local into the (x, y) coodinate system.
/** Let n be the normal and tau the tangential. Then this computes
 *
 *  u_rot = (rho, u[1] * n + u[2] * tau, E)
 *
 */
EulerState globalCoordinates(const EulerState &u,
                             const Eigen::Vector2d &normal) {
    EulerState u_rot;

    u_rot[0] = u[0];

    u_rot[1] = u[1] * normal[0] - u[2] * normal[1];
    u_rot[2] = u[1] * normal[1] + u[2] * normal[0];

    u_rot[3] = u[3];

    return u_rot;
}

/// Compute the conserved variables given the primitives.
EulerState conservedVars(const EulerState &w) {
    double rho = w[0];
    double vx = w[1];
    double vy = w[2];
    double p = w[3];

    EulerState u;
    u[0] = rho;
    u[1] = rho * vx;
    u[2] = rho * vy;
    u[3] = energyFromPVars(rho, vx, vy, p);

    return u;
}

/// Compute the primitive variables given the conserved ones.
EulerState primitiveVars(const EulerState &u) {
    double rho = u[0];
    double mx = u[1];
    double my = u[2];

    EulerState w;
    w[0] = rho;
    w[1] = mx / rho;
    w[2] = my / rho;
    w[3] = pressure(u);

    return w;
}

}
