#pragma once

/// HLLC numerical flux with Einfeldt-Batten wavespeeds.
/** Reference: Batten, Wavespeed Estimates for the HLLC Riemann Solver, 1997
 *  @param euler
 *  @param uL    conserved variables 'inside' of the cell.
 *  @param uR    conserved variables 'outside' of the cell.
 */
EulerState hllc(const EulerState &uL, const EulerState &uR) {
  // implement the HLLC approximate Riemann solver.
  // tip, compute the three wave speeds sL, s_star & sR in
  // a separate function to keep things readable.

  // convert conservative to primitive
  const auto uL_prim = euler::primitiveVars(uL);
  const auto uR_prim = euler::primitiveVars(uR);

  // extract data, i.e. rho, v, w, p, E
  // v -> x-component of speed
  // z -> y-component of speed
  const double rhoL = uL_prim(0);
  const double rhoR = uR_prim(0);
  const double vL = uL_prim(1);
  const double vR = uR_prim(1);
  const double zL = uL_prim(2);
  const double zR = uR_prim(2);
  const double pL = uL_prim(3);
  const double pR = uR_prim(3);
  const double EL = uL(3);
  const double ER = uR(3);

  /* speeds */
  const double sL =
      euler::minEigenValue(uL, uR); // TODO: Implement Einfeldt-Batten sL
  const double sR =
      euler::maxEigenValue(uL, uR); // TODO: Implement Einfeldt-Batten sR
  const double s_star =
      (pR - pL + rhoL * vL * (sL - vL) - rhoR * vR * (sR - vR)) /
      (rhoL * (sL - vL) - rhoR * (sR - vR));

  /* uL_star state */
  EulerState uL_star;
  uL_star(0) = 1.0;
  uL_star(1) = s_star;
  uL_star(2) = zL;
  uL_star(3) =
      (EL / rhoL) + (s_star - vL) * (s_star + (pL / (rhoL * (sL - vL))));
  uL_star = rhoL * ((sL - vL) / (sL - s_star)) * uL_star;

  /* uR_star state */
  EulerState uR_star;
  uR_star(0) = 1.0;
  uR_star(1) = s_star;
  uR_star(2) = zR;
  uR_star(3) =
      (ER / rhoR) + (s_star - vR) * (s_star + (pR / (rhoR * (sR - vR))));
  uR_star = rhoR * ((sR - vR) / (sR - s_star)) * uR_star;

  /* HLLC Fluxes */
  const auto fL = euler::flux(uL);
  const auto fR = euler::flux(uR);
  const auto fL_star = fL + sL * (uL_star - uL);
  const auto fR_star = fR + sR * (uR_star - uR);

  if (sL > 0)
    return fL;
  else if (sL <= 0 && s_star > 0)
    return fL_star;
  else if (s_star <= 0 && sR > 0)
    return fR_star;
  else /* if ( sR <= 0) */
    return fR;
}
