#pragma once

double sign(double a) { return std::copysign(1.0, a); }

double minmod(double a, double b) {
  return 0.5 * (sign(a) + sign(b)) * std::min(std::abs(a), std::abs(b));
}

double maxmod(double a, double b) {
  return 0.5 * (sign(a) + sign(b)) * std::max(std::abs(a), std::abs(b));
}

double minmod(double a, double b, double c) { return minmod(a, minmod(b, c)); }

// Uses minmod unless having explicitly defined superbee or mc
double slope_limiter(double sL, double sR) {
#if defined(ANCSE_SUPERBEE_LIMITER)
  return maxmod(minmod(2.0 * sL, sR), minmod(sL, 2.0 * sR));
#elif defined(ANCSE_MC_LIMITER)
  return minmod(2.0 * sL, 0.5 * (sL + sR), 2.0 * sR);
#else
  return minmod(sL, sR);
#endif
}
