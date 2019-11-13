#include "solve_euler.hpp"
#include "writer.hpp"

int main(int, char **) {

    auto mesh = load_naca_mesh(ANCSE_DATA_PATH "/naca.mesh");

    int n_cells = mesh.getNumberOfTriangles();
    Eigen::MatrixXd u0(n_cells, 4);

    for (int i = 0; i < n_cells; ++i) {
        double T = 1.0;
        double p = 1.0;
        double cs = std::sqrt(GAMMA * T);
        double M = 0.85;    // Mach number
        double alpha = 1.0; // angle of attack

        double rho = p / T;
        double vx = M * cs * std::cos(alpha * M_PI / 180.0);
        double vy = M * cs * std::sin(alpha * M_PI / 180.0);

        u0(i, 0) = rho;
        u0(i, 1) = rho * vx;
        u0(i, 2) = rho * vy;
        u0(i, 3) = euler::energyFromPVars(rho, vx, vy, p);
    }

    auto U = solveEuler(u0, mesh, 100.0);

    writeMatrixToFile("naca_u.txt", U);
    mesh.writeToFile("naca");
}
