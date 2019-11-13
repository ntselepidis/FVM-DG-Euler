#ifndef HYPSYS1D_TIME_LOOP_HPP
#define HYPSYS1D_TIME_LOOP_HPP

#include <Eigen/Dense>
#include <fmt/format.h>
#include <fstream>

#include <ancse/cfl_condition.hpp>
#include <ancse/runge_kutta.hpp>
#include <ancse/simulation_time.hpp>
#include <ancse/snapshot_writer.hpp>

/// Advances the solution of the PDE iteratively until the final time.
class TimeLoop {
  public:
    TimeLoop(std::shared_ptr<SimulationTime> simulation_time,
             std::shared_ptr<TimeIntegrator> time_integrator,
             std::shared_ptr<CFLCondition> cfl_condition,
             std::shared_ptr<SnapshotWriter> snapshot_writer)
        : simulation_time(std::move(simulation_time)),
          time_integrator(std::move(time_integrator)),
          cfl_condition(std::move(cfl_condition)),
          snapshot_writer(std::move(snapshot_writer)) {}

    virtual ~TimeLoop() {}

    virtual void operator()(Eigen::MatrixXd u0) const;

    void write_snapshot(const Eigen::MatrixXd &u) const;

  protected:
    mutable std::shared_ptr<SimulationTime> simulation_time;
    std::shared_ptr<TimeIntegrator> time_integrator;
    std::shared_ptr<CFLCondition> cfl_condition;
    std::shared_ptr<SnapshotWriter> snapshot_writer;
};

#endif // HYPSYS1D_TIME_LOOP_HPP
