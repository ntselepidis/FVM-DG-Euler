// Putting the header as the very first line of code, ensures that your
// header is self-contained. It's impossible for include-ordering to
// break your code.
#include <ancse/snapshot_writer.hpp>

#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

JSONSnapshotWriter <FVM>
:: JSONSnapshotWriter(const Grid &grid_,
                      std::shared_ptr<Model> &model_,
                      std::shared_ptr<SimulationTime> simulation_time_,
                      std::string basedir_,
                      std::string basename_)
    : grid(grid_),
      model(model_),
      simulation_time(std::move(simulation_time_)),
      basedir(std::move(basedir_)),
      basename(std::move(basename_))
{
    fs::create_directories(basedir);
}

void JSONSnapshotWriter <FVM>
:: operator()(const Eigen::MatrixXd &u) const
{
    std::string filename = basedir + "/" + basename
                         + fmt::format("{:04d}.json", k_output);
    k_output += 1;

    auto json = nlohmann::json{};
    unsigned int n_cells = static_cast<unsigned int> (grid.n_cells);

    if ((model->get_name()).compare("burgers") == 0)
    {
        for (unsigned int i = 0; i < n_cells; ++i) {
            json["data"][i] = u(0,i);
        }
    }
    else if ((model->get_name()).compare("euler") == 0)
    {
        for (unsigned int i = 0; i < n_cells; ++i) 
        {
            auto u_prim = model->cons_to_prim(u.col(i));

            json["density"][i]  = u_prim(0);
            json["velocity"][i] = u_prim(1);
            json["pressure"][i] = u_prim(2);
        }
    }

    for (unsigned int i = 0; i < n_cells; ++i) {
        json["cell_centers"][i] = cell_center(grid, static_cast<int>(i));
    }
    json["time"] = simulation_time->t;

    auto file = std::ofstream(filename);
    assert(file.good());  // always check that you can write to the file.

    file << json.dump(2);
}


JSONSnapshotWriter <DG>
:: JSONSnapshotWriter(const Grid &grid_,
                      std::shared_ptr<Model> &model_,
                      const DGHandler &dg_handler_,
                      std::shared_ptr<SimulationTime> simulation_time_,
                      std::string basedir_,
                      std::string basename_)
    : grid(grid_),
      model(model_),
      dg_handler(dg_handler_),
      simulation_time(std::move(simulation_time_)),
      basedir(std::move(basedir_)),
      basename(std::move(basename_))
{
    fs::create_directories(basedir);
}

void JSONSnapshotWriter <DG>
:: operator()(const Eigen::MatrixXd &u) const
{
    std::string filename = basedir + "/" + basename
                         + fmt::format("{:04d}.json", k_output);
    k_output += 1;

    auto json = nlohmann::json{};
    unsigned int n_cells = static_cast<unsigned int> (grid.n_cells);


    if ((model->get_name()).compare("burgers") == 0)
    {
        auto uSol = dg_handler.build_cell_avg(u);
        for (unsigned int i = 0; i < n_cells; ++i) {
            json["data"][i] = uSol(0,i);
        }
    }
    else if ((model->get_name()).compare("euler") == 0)
    {
        auto uSol = dg_handler.build_cell_avg(u);
        for (unsigned int i = 0; i < n_cells; ++i) 
        {
            auto u_prim = model->cons_to_prim(uSol.col(i));
        
            json["density"][i]  = u_prim(0);
            json["velocity"][i] = u_prim(1);
            json["pressure"][i] = u_prim(2);
        }
    }

    for (unsigned int i = 0; i < n_cells; ++i) {
        json["cell_centers"][i] = cell_center(grid, static_cast<int>(i));
    }
    json["time"] = simulation_time->t;

    auto file = std::ofstream(filename);
    assert(file.good());  // always check that you can write to the file.

    file << json.dump(2);
}
