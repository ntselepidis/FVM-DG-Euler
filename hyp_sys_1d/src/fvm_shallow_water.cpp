#include <Eigen/Dense>

#include <ancse/config.hpp>
#include <ancse/cfl_condition.hpp>
#include <ancse/fvm_rate_of_change.hpp>
#include <ancse/snapshot_writer.hpp>
#include <ancse/time_loop.hpp>

static int n_vars = 2;

template<class F>
Eigen::MatrixXd ic(const F &f, const Grid &grid) {
    Eigen::MatrixXd u0(n_vars, grid.n_cells);
    for(int i = 0; i < grid.n_cells; ++i) {
        u0.col(i) = f(cell_center(grid, i));
    }

    return u0;
}

TimeLoop make_fvm(const nlohmann::json &config,
                  const Grid &grid)
{
    double t_end = config["t_end"];
    double cfl_number = config["cfl_number"];

    auto n_ghost = grid.n_ghost;
    auto n_cells = grid.n_cells;

    std::shared_ptr<Model> model = make_model(config);
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
    auto snapshot_writer = std::make_shared<JSONSnapshotWriter<FVM>>
            (grid, model, simulation_time,
             std::string(config["output_dir"]),
             std::string(config["output_file"]));

    return TimeLoop(simulation_time, time_integrator,
                    cfl_condition, snapshot_writer);
}

void dam_break_test(const nlohmann::json &config)
{
    auto fn = [](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0) { // left state
            u(0) = 3;
            u(1) = 0;
        } else {        // right state
            u(0) = 1;
            u(1) = 0;
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({-1, 1}, n_cells, n_ghost);
    auto u0 = ic(fn, grid);
    auto fvm = make_fvm(config, grid);
    fvm(u0);
}

void vacuum_test(const nlohmann::json &config)
{
    auto fn = [](double x) {
        Eigen::VectorXd u(n_vars);
        if (x <= 0) { // left state
            u(0) = 1;
            u(1) = -1;
        } else {        // right state
            u(0) = 1;
            u(1) = +1;
        }
        return u;
    };

    int n_ghost = config["n_ghost"];
    int n_cells = int(config["n_interior_cells"]) + n_ghost * 2;

    auto grid = Grid({-1, 1}, n_cells, n_ghost);
    auto u0 = ic(fn, grid);
    auto fvm = make_fvm(config, grid);
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
    if (ic_key == "dam_break") {
        dam_break_test(config);
    } else if (ic_key == "vacuum") {
        vacuum_test(config);
    }

    return 0;
}
