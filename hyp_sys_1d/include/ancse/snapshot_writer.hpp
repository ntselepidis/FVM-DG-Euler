#ifndef HYPSYS1D_SNAPSHOT_WRITER_HPP
#define HYPSYS1D_SNAPSHOT_WRITER_HPP

#include <Eigen/Dense>
#include <memory>

#include <ancse/includes.hpp>
#include <ancse/grid.hpp>
#include <ancse/simulation_time.hpp>
#include <ancse/model.hpp>
#include <ancse/dg_handler.hpp>

/// Interface for outputting snapshots.
class SnapshotWriter {
  public:
    virtual ~SnapshotWriter() = default;

    virtual void operator()(const Eigen::MatrixXd &u) const = 0;
};

template <int DiscretisationType>
class JSONSnapshotWriter;

template<>
class JSONSnapshotWriter <FVM> : public SnapshotWriter {
  public:
    JSONSnapshotWriter(const Grid &grid,
                       std::shared_ptr<Model> &model,
                       std::shared_ptr<SimulationTime> simulation_time,
                       std::string basedir,
                       std::string basename);

    virtual void operator()(const Eigen::MatrixXd &u) const override;

  private:
    Grid grid;
    std::shared_ptr<Model> model;
    std::shared_ptr<SimulationTime> simulation_time;

    std::string basedir;
    std::string basename;
    mutable int k_output = 0;
};

template<>
class JSONSnapshotWriter <DG> : public SnapshotWriter {
  public:
    JSONSnapshotWriter(const Grid &grid,
                       std::shared_ptr<Model> &model,
                       const DGHandler &dg_handler,
                       std::shared_ptr<SimulationTime> simulation_time,
                       std::string basedir,
                       std::string basename);

    virtual void operator()(const Eigen::MatrixXd &u) const override;

  private:
    Grid grid;
    std::shared_ptr<Model> model;
    DGHandler dg_handler;
    std::shared_ptr<SimulationTime> simulation_time;

    std::string basedir;
    std::string basename;
    mutable int k_output = 0;
};

#endif // HYPSYS1D_SNAPSHOT_WRITER_HPP
