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
    return EulerState{};
}
