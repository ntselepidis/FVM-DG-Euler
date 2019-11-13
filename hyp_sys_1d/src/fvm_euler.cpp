#include <Eigen/Dense>

#include <ancse/config.hpp>
#include <ancse/cfl_condition.hpp>
#include <ancse/fvm_rate_of_change.hpp>
#include <ancse/snapshot_writer.hpp>
#include <ancse/time_loop.hpp>

static int n_vars = 3;

template<class F>
Eigen::MatrixXd ic(const F &f, const Grid &grid) {
    Eigen::MatrixXd u0(n_vars, grid.n_cells);
    for(int i = 0; i < grid.n_cells; ++i) {
        u0.col(i) = f(cell_center(grid, i));
    }

    return u0;
}

TimeLoop make_fvm(const nlohmann::json &config,
                  const Grid &grid,
                  std::shared_ptr<Model> &model)
{
    double t_end = config["t_end"];
    double cfl_number = config["cfl_number"];

    auto n_ghost = grid.n_ghost;
    auto n_cells = grid.n_cells;

    auto n_vars = model->get_nvars();

    auto simulation_time = std::make_shared<SimulationTime>(t_end);
    auto fvm_rate_of_change
            = make_fvm_rate_of_change(config, grid, model,
                                      simulation_time);
    auto boundary_condition
            = make_boundary_condition(n_ghost,
                                      config["boundary_condition"]);
    auto time_integrator = make_runge_kutta(config,
                                            fvm_rate_of_change,
                                            boundary_condition,
                                            n_vars, n_cells);
    auto cfl_condition = make_cfl_condition(grid, model, cfl_number);
    auto snapshot_writer = std::make_shared< JSONSnapshotWriter<FVM> >
            (grid, model, simulation_time,
             std::string(config["output_dir"]),
             std::string(config["output_file"]));

    return TimeLoop(simulation_time, time_integrator,
                    cfl_condition, snapshot_writer);
}

void sod_shock_tube_test(const nlohmann::json &config)
{
    double gamma = 7./5.;
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto model_euler = dynamic_cast<Euler*>(model.get());
    model_euler->set_gamma(gamma);

    auto fn = [gamma](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0.5) { // left state
            u(0) = 1;
            u(1) = 0;
            u(2) = 1/(gamma-1);
        } else {        // right state
            u(0) = 0.125;
            u(1) = 0;
            u(2) = 0.1/(gamma-1);
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({0.0, 1.0}, n_cells, n_ghost);
    auto u0 = ic(fn, grid);

    auto fvm = make_fvm(config, grid, model);
    fvm(u0);
}

void vacuum_test(const nlohmann::json &config)
{
    double gamma = 7./5.;
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto model_euler = dynamic_cast<Euler*>(model.get());
    model_euler->set_gamma(gamma);

    auto fn = [](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0.5) { // left state
            u(0) = 1;
            u(1) = -2;
            u(2) = 4.5;
        } else {        // right state
            u(0) = 1;
            u(1) = 2;
            u(2) = 4.5;
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({0.0, 1.0}, n_cells, n_ghost);
    auto u0 = ic(fn, grid);

    auto fvm = make_fvm(config, grid, model);
    fvm(u0);
}

void lax_shock_tube_test(const nlohmann::json &config)
{
    double gamma = 7./5.;
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto model_euler = dynamic_cast<Euler*>(model.get());
    model_euler->set_gamma(gamma);

    auto fn = [](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0.5) { // left state
            u(0) = 0.445;
            u(1) = 0.311;
            u(2) = 8.928;
        } else {        // right state
            u(0) = 0.5;
            u(1) = 0;
            u(2) = 1.4275;
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({0.0, 1.0}, n_cells, n_ghost);
    auto u0 = ic(fn, grid);

    auto fvm = make_fvm(config, grid, model);
    fvm(u0);
}

void smooth_test(const nlohmann::json &config)
{
    double gamma = 7./5.;
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto model_euler = dynamic_cast<Euler*>(model.get());
    model_euler->set_gamma(gamma);

    double a = 0.1;
    double sigma = 0.1;
    auto fn = [gamma, a, sigma](double x) {
        Eigen::VectorXd u(n_vars);
        u(0) = 1;
        u(1) = 0;
        u(2) = (1 + a*exp(-x*x/(sigma*sigma)))/(gamma-1);
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({-1.0, +1.0}, n_cells, n_ghost);
    auto u0 = ic(fn, grid);

    auto fvm = make_fvm(config, grid, model);
    fvm(u0);
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
    if (ic_key == "sod_shock_tube") {
        sod_shock_tube_test(config);
    } else if (ic_key == "vacuum") {
    vacuum_test(config);
    } else if (ic_key == "lax_shock_tube") {
        lax_shock_tube_test(config);
    } else if (ic_key == "smooth") {
        smooth_test(config);
    }

    return 0;
}
