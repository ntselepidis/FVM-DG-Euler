#include <ancse/dg_handler.hpp>


/// build solution from DG coefficients and the basis
/// pre-evaluated at a certain point
Eigen::VectorXd DGHandler
:: build_sol(const Eigen::VectorXd& u,
             const Eigen::VectorXd& basis) const
{
    Eigen::VectorXd uSol(n_vars);

    for (int i = 0; i < n_vars; i++) {
        uSol(i) = u.segment(i*n_coeff, n_coeff).dot(basis);
    }

    return uSol;
}

/// build solution from DG coefficients at a given reference point
Eigen::VectorXd DGHandler 
:: build_sol(const Eigen::VectorXd& u,
             double xi) const
{
    Eigen::VectorXd uSol(n_vars);

    uSol = build_sol(u, poly_basis(xi));

    return uSol;
}

/// build cell average
Eigen::MatrixXd DGHandler
:: build_cell_avg (const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();
    Eigen::MatrixXd u0 (n_vars, n_cells);

    // Cell avg is eq to 0-th ord coeffs of each component (scaled)
    for (int j = 0; j < n_cells; j++) {
        for (int i = 0; i < n_vars; i++) {
            u0(i,j) = u(i*n_coeff,j);
        }
    }

    // scaling factor ( 1.0 / sqrt(h) )
    const double scal = poly_basis(0.0)(0);

    return scal*u0;
}

/// build split solution uSol_m = u0 + um, uSol_p = u0 - up
/// from DG coefficients
std::tuple <Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
DGHandler :: build_split_sol(const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();
    if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    
    auto u0 = build_cell_avg(u);
    Eigen::MatrixXd um(n_vars, n_cells);
    Eigen::MatrixXd up(n_vars, n_cells);
    auto phim = poly_basis(1.0);
    auto phip = poly_basis(0.0);

    for (int j = 0; j < n_cells; j++)
    {
        for (int i = 0; i < n_vars; i++)
        {
            um(i,j) = u.col(j).segment(i*n_coeff+1, n_coeff-1)
                            .dot(phim.segment(1, n_coeff-1));

            up(i,j) =-u.col(j).segment(i*n_coeff+1, n_coeff-1)
                            .dot(phip.segment(1, n_coeff-1));
        }
    }

    return {std::move(u0), std::move(um), std::move(up)};
}

/// build DG coefficients from uSol_m = u0 + um, uSol_p = u0 - up
void DGHandler :: compute_limit_coeffs (Eigen::MatrixXd &u,
                                        Eigen::MatrixXd &um,
                                        Eigen::MatrixXd &up) const
{
    if (n_coeff == 1) {
        return;
    }
    else if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    



}
