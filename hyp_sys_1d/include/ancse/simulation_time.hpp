#ifndef HYPSYS1D_SIMULATION_TIME_HPP
#define HYPSYS1D_SIMULATION_TIME_HPP

/// Time, step-size and number of time-steps of the simulation.
/** In FVM codes there is occasionally the need to access to either the
 *  time-step 'dt' or iteration count 'k' in unpredicted places.
 */
struct SimulationTime {
    explicit SimulationTime(double t_end);
    SimulationTime(double t, double dt, double t_end, int k);

    double t;     ///< current time of the simulation
    double dt;    ///< size of the time-step
    double t_end; ///< simulation runs until `t_end`.
    int k;        ///< number of time-steps performed

    void advance();
};

bool is_finished(const SimulationTime &simulation_time);

#endif // HYPSYS1D_SIMULATION_TIME_HPP
