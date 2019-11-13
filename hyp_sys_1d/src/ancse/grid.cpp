#include <ancse/grid.hpp>

Grid::Grid(std::pair<double, double> domain, int n_cells, int n_ghost)
    : domain(domain), n_cells(n_cells), n_ghost(n_ghost) {
    auto [a, b] = domain;
    dx = (b - a) / (n_cells - 2 * n_ghost);
}

double cell_center(const Grid &grid, int i) {
    auto [a, _] = grid.domain;
    return a + grid.dx * (i + 0.5 - grid.n_ghost);
}

double cell_point(const Grid &grid, int i, double xi) {
    auto [a, _] = grid.domain;
    return a + grid.dx * (i + xi - grid.n_ghost);
}
