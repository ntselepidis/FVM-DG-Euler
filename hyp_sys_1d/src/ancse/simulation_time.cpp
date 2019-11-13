#include <ancse/simulation_time.hpp>

SimulationTime::SimulationTime(double t_end)
    : t(0.0), dt(-1.0), t_end(t_end), k(0) {}

SimulationTime::SimulationTime(double t, double dt, double t_end, int k)
    : t(t), dt(dt), t_end(t_end), k(k) {}

void SimulationTime::advance() {
    t = t + dt;
    k += 1;
}
bool is_finished(const SimulationTime &simulation_time) {
    return simulation_time.t
           > simulation_time.t_end - 1e-12 * simulation_time.t_end;
}
