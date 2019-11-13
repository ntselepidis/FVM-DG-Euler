#pragma once

#include "cfl_condition.hpp"
#include "mesh.hpp"
#include "numerical_flux.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>
#include <fmt/format.h>

void print_progress(
    int n_steps,
    double T,
    double t,
    double dt,
    const std::chrono::high_resolution_clock::time_point &time_start);

Eigen::MatrixXd solveEuler(Eigen::MatrixXd U0, const Mesh &mesh, double T) {

    int n_cells = mesh.getNumberOfTriangles();

    Eigen::MatrixXd U_tmp(n_cells, 4);
    Eigen::MatrixXd dUdt(n_cells, 4);

    double t = 0;

    auto cfl_condition = CFLCondition(mesh);
    auto computeNetFlux = FluxRateOfChange(n_cells);

    int n_steps = 0;

    auto time_start = std::chrono::high_resolution_clock::now();

    std::cout << std::endl;
    while (t < T) {
        double dt = cfl_condition(U0);

        // Compute u^(*)
        computeNetFlux(dUdt, U0, mesh);
#pragma omp parallel for
        for (int i = 0; i < n_cells; ++i) {
            U_tmp.row(i) = U0.row(i) + dt * dUdt.row(i);
        }

        // Compute u^(**) and u^{n+1} in one step. We can update U0 directly.
        computeNetFlux(dUdt, U_tmp, mesh);
#pragma omp parallel for
        for (int i = 0; i < n_cells; ++i) {
            U0.row(i) = 0.5 * (U0.row(i) + U_tmp.row(i) + dt * dUdt.row(i));
        }

        t += dt;
        n_steps++;

        print_progress(n_steps, T, t, dt, time_start);
    }
    std::cout << "\n";

    return U0;
}

void print_progress(
    int n_steps,
    double T,
    double t,
    double dt,
    const std::chrono::high_resolution_clock::time_point &time_start) {

    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(
        now - time_start);
    double elapsed = duration.count();
    double eta = elapsed * (T - t) / t;

    std::cout << fmt::format("\r"
                             "n_steps = {: 3d}, t = {:.3e}, dt = {:.2e}, "
                             "elased = {:6.1f}s, eta = {:6.1f}s",
                             n_steps,
                             t,
                             dt,
                             elapsed,
                             eta);
    std::cout.flush();
}
