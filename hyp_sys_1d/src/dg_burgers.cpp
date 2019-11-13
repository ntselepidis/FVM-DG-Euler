#include <Eigen/Dense>
#include <iostream>

#include <ancse/includes.hpp>
#include <ancse/config.hpp>
#include <ancse/polynomial_basis.hpp>
#include <ancse/dg_handler.hpp>
#include <ancse/dg_limiting.hpp>
#include <ancse/cfl_condition.hpp>
#include <ancse/dg_rate_of_change.hpp>
#include <ancse/snapshot_writer.hpp>
#include <ancse/time_loop.hpp>

static int n_vars = 1;

template<class F>
Eigen::MatrixXd ic(const F &f,
                   const Grid &grid,
                   const PolynomialBasis &poly_basis,
                   const DGHandler &dg_handler)
{
    int n_coeff = 1 + poly_basis.get_degree();
    Eigen::MatrixXd u0 = Eigen::MatrixXd::Zero (n_vars*n_coeff,
                                                grid.n_cells);

    auto [quad_points, quad_weights] = dg_handler.get_quadrature();
    int n_quad = static_cast<int>(quad_points.size());

    /// eval basis and its derivate for all quadrature points
    Eigen::MatrixXd basis(n_coeff, n_quad);
    for (int k = 0; k < n_quad; k++) {
        basis.col(k) = poly_basis(quad_points(k));
    }

    // L2-projection
    for (int j = 0; j < grid.n_cells; ++j)
    {
        for (int k = 0; k < quad_points.size(); k++)
        {
            auto fVal = f(cell_point(grid, j, quad_points(k)));
            for (int i = 0; i < n_vars; i++) {
                u0.col(j).segment(i*n_coeff, n_coeff)
                        += quad_weights(k)*fVal(i)*basis.col(k);
            }
        }
    }

    return u0*grid.dx;
}

TimeLoop make_dg(const nlohmann::json &config,
                 const Grid &grid,
                 const PolynomialBasis &poly_basis,
                 const DGHandler &dg_handler,
                 std::shared_ptr<Model> &model)
{
    double t_end = config["t_end"];
    double cfl_number = config["cfl_number"];

    auto n_ghost = grid.n_ghost;
    auto n_cells = grid.n_cells;

    auto n_vars = model->get_nvars();
    int n_coeff = 1 + poly_basis.get_degree();

    auto simulation_time = std::make_shared<SimulationTime>(t_end);
    auto boundary_condition
            = make_boundary_condition(n_ghost,
                                      config["boundary_condition"]);
    auto dg_rate_of_change
            = make_dg_rate_of_change(config, grid, model,
                                     poly_basis, dg_handler,
                                     simulation_time);
    auto dg_limiting
            = make_dg_limiting(config, grid, dg_handler);
    auto time_integrator = make_runge_kutta(config,
                                            dg_rate_of_change,
                                            boundary_condition,
                                            dg_limiting,
                                            n_vars*n_coeff, n_cells);
    auto cfl_condition
            = make_cfl_condition(grid, model, dg_handler, cfl_number);
    auto snapshot_writer = std::make_shared<JSONSnapshotWriter <DG>>
            (grid, model, dg_handler, simulation_time,
             std::string(config["output_dir"]),
             std::string(config["output_file"]));

    return TimeLoop(simulation_time, time_integrator,
                    cfl_condition, snapshot_writer);
}

void shock_test(const nlohmann::json &config)
{
    int deg = int(config["degree"]);

    std::shared_ptr<Model> model = std::make_shared<Burgers>();

    auto fn = [](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0.5) { // left state
            u(0) = 1;
        } else {        // right state
            u(0) = 0;
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({0.0, 1.0}, n_cells, n_ghost);
    auto poly_basis = PolynomialBasis(deg, 1./sqrt(grid.dx));
    auto dg_handler = DGHandler(model, poly_basis);
    auto u0 = ic(fn, grid, poly_basis, dg_handler);

    auto dg = make_dg(config, grid, poly_basis, dg_handler, model);
    dg(u0);
}

void rarefaction_test(const nlohmann::json &config)
{
    int deg = int(config["degree"]);

    std::shared_ptr<Model> model = std::make_shared<Burgers>();

    auto fn = [](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0.5) { // left state
            u(0) = 0;
        } else {        // right state
            u(0) = 1;
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({0.0, 1.0}, n_cells, n_ghost);
    auto poly_basis = PolynomialBasis(deg, 1./sqrt(grid.dx));
    auto dg_handler = DGHandler(model, poly_basis);
    auto u0 = ic(fn, grid, poly_basis, dg_handler);

    auto dg = make_dg(config, grid, poly_basis, dg_handler, model);
    dg(u0);
}

int main(int argc, char* const argv[])
{
    nlohmann::json config;
    std::string fileName;
    if (argc == 2) {
        fileName = argv[1];
    } else {
        fileName = "../config.json";
    }
    config = get_config (fileName);

    std::string ic_key = config["initial_conditions"];
    if (ic_key == "shock") {
        shock_test(config);
    } else if (ic_key == "rarefaction") {
        rarefaction_test(config);
    }

    return 0;
}
