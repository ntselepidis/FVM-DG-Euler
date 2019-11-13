#include<gtest/gtest.h>

#include <ancse/boundary_condition.hpp>


TEST(TestBoundaryCondition, Periodic) {
    int n_cells = 10;
    int n_ghost = 2;
    
    int n_vars = 2;
    Eigen::MatrixXd u(n_vars, n_cells);
    for(int j = 0; j < n_cells; ++j) {
        for(int i = 0; i < n_vars; ++i) {
            u(i,j) = i + j*n_vars;
        }
    }

    auto bc = PeriodicBC(n_ghost);
    bc(u);

    ASSERT_DOUBLE_EQ(u(0,0), 12.0);
    ASSERT_DOUBLE_EQ(u(1,0), 13.0);
    ASSERT_DOUBLE_EQ(u(0,1), 14.0);
    ASSERT_DOUBLE_EQ(u(1,1), 15.0);

    for(int j = n_ghost; j < n_cells-n_ghost; ++j) {
        for(int i = 0; i < n_vars; ++i) {
            ASSERT_DOUBLE_EQ(u(i,j), i + j*n_vars) << "Failed on cell = " << j << " , at var = " << i;
        }
    }
    
    ASSERT_DOUBLE_EQ(u(0,8), 4.0);
    ASSERT_DOUBLE_EQ(u(1,8), 5.0);
    ASSERT_DOUBLE_EQ(u(0,9), 6.0);
    ASSERT_DOUBLE_EQ(u(1,9), 7.0);
}

TEST(TestBoundaryCondition, Outflow) {
    int n_cells = 10;
    int n_ghost = 2;

    int n_vars = 2;
    Eigen::MatrixXd u(n_vars, n_cells);
    for(int j = 0; j < n_cells; ++j) {
        for(int i = 0; i < n_vars; ++i) {
            u(i,j) = i + j*n_vars;
        }
    }

    auto bc = OutflowBC(n_ghost);
    bc(u);
    
    ASSERT_DOUBLE_EQ(u(0,0), 4.0);
    ASSERT_DOUBLE_EQ(u(1,0), 5.0);
    ASSERT_DOUBLE_EQ(u(0,1), 4.0);
    ASSERT_DOUBLE_EQ(u(1,1), 5.0);

    for(int j = n_ghost; j < n_cells-n_ghost; ++j) {
        for(int i = 0; i < n_vars; ++i) {
            ASSERT_DOUBLE_EQ(u(i,j), i + j*n_vars) << "Failed on cell = " << j << " , at var = " << i;
        }
    }
    
    ASSERT_DOUBLE_EQ(u(0,8), 14.0);
    ASSERT_DOUBLE_EQ(u(1,8), 15.0);
    ASSERT_DOUBLE_EQ(u(0,9), 14.0);
    ASSERT_DOUBLE_EQ(u(1,9), 15.0);
}
