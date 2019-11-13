#include <gtest/gtest.h>

#include <ancse/grid.hpp>

TEST(Grid, cell_center) {
    double a = 2.0, b = 4.0;
    double dx = 2.0 / 10.0;

    int n_cells = 14;
    int n_ghost = 2;

    auto grid = Grid({a, b}, n_cells, n_ghost);

    // The cell-center of the first interior cell is
    //    a + 0.5 *dx
    ASSERT_DOUBLE_EQ(cell_center(grid, n_ghost), a + 0.5*dx);
}
