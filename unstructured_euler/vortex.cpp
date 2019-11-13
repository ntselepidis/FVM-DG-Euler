#include "solve_euler.hpp"
#include "writer.hpp"
#include <fmt/format.h>

void run_vortex_simulation(const Mesh &mesh, const std::string &basename) {
    auto n_cells = mesh.getNumberOfTriangles();

    Eigen::MatrixXd u0(n_cells, 4);

    for (int i = 0; i < n_cells; ++i) {
        auto x = mesh.getCellCenter(i);
        auto xc = Eigen::Vector2d{};
        xc[0] = 0.5;
        xc[1] = 0.5;

        double r = (x - xc).norm();
        double rc = 0.1;
        double eps = 5.0 / (2.0 * M_PI);
        double alpha = 0.5;
        double tau = r / rc;

        double phi_r = eps * std::exp(alpha * (1.0 - tau * tau));
        double vx = -(x[1] - xc[1]) * phi_r;
        double vy = +(x[0] - xc[0]) * phi_r;
        double theta = 1.0 - (GAMMA - 1.0) / (2.0 * GAMMA) * phi_r * phi_r;

        // theta = p / rho
        // p = rho ** gamma
        // =>
        // theta = rho ** (gamma - 1)
        // rho = theta ** (-gamma + 1)
        double rho = std::pow(theta, 1.0 / (GAMMA - 1.0));

        //        double p = theta * rho;
        double p = 0.1 + std::exp(-(r * r / 0.1));

        u0(i, 0) = rho;
        u0(i, 1) = rho * vx;
        u0(i, 2) = rho * vy;
        u0(i, 3) = p / (GAMMA - 1) + 0.5 * rho * (vx * vx + vy * vy);
    }

    auto u1 = solveEuler(u0, mesh, 0.02);

    writeMatrixToFile(basename + "_u.txt", u1);
    mesh.writeToFile(basename);
}

int main(int, char **) {

    std::vector<std::pair<std::string, std::string>> experiments;
    for (int k = 1; k <= 6; ++k) {
        auto filename = ANCSE_DATA_PATH + fmt::format("/square_{}.mesh", k);
        auto experiment = fmt::format("vortex_{}", k);
        experiments.emplace_back(filename, experiment);
    }

    for (const auto &[mesh_name, basename] : experiments) {
        run_vortex_simulation(load_square_mesh(mesh_name), basename);
    }
}
