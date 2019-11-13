#include <ancse/time_loop.hpp>

#include <iostream>


void TimeLoop :: operator()(Eigen::MatrixXd u0) const
{
    Eigen::MatrixXd u1(u0.rows(), u0.cols());

    simulation_time->dt = (*cfl_condition)(u0);

    while (!is_finished(*simulation_time)) {
        (*time_integrator)(u1, u0, simulation_time->dt);

        simulation_time->advance();
        
        double dt_cfl = (*cfl_condition)(u1);
        simulation_time->dt = dt_cfl;
        if (simulation_time->t + dt_cfl > simulation_time->t_end) {
            simulation_time->dt = simulation_time->t_end - simulation_time->t;
        }

        //write_snapshot(u1);

        u0.swap(u1);

        std::cout << fmt::format("{: 4d}: {} {}\n",
                                 simulation_time->k,
                                 simulation_time->t,
                                 simulation_time->dt);
    }
    write_snapshot(u0);
}

void TimeLoop :: write_snapshot(const Eigen::MatrixXd &u) const {
    (*snapshot_writer)(u);
}
