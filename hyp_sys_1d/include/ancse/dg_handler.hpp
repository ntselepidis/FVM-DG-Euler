#ifndef HYPSYS1D_DG_HANDLER_HPP
#define HYPSYS1D_DG_HANDLER_HPP

#include <memory>
#include <Eigen/Dense>
#include <fmt/format.h>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/polynomial_basis.hpp>


class DGHandler {
  public:

    DGHandler (const std::shared_ptr<Model> &model,
               const PolynomialBasis &poly_basis)
        : poly_basis (poly_basis)
    {
        n_vars = model->get_nvars();
        n_coeff = 1 + poly_basis.get_degree();
        
        n_quad = 3;
        /// Gauss-Legendre quadrature, degree 5
        quad_points.resize(n_quad);
        quad_points << (1 - sqrt(3./5.))/2., 0.5, (1 + sqrt(3./5.))/2.;

        quad_weights.resize(n_quad);
        quad_weights << 5./18., 4./9., 5./18.;
    }

    /// set quadrature
    void set_quadrature (const Eigen::VectorXd &quad_points_,
                         const Eigen::VectorXd &quad_weights_) const
    {
        quad_points = quad_points_;
        quad_weights = quad_weights_;
    }
    
    /// get quadrature
    std::pair <Eigen::VectorXd, Eigen::VectorXd> get_quadrature () const 
    {
        return {quad_points, quad_weights};
    }

    /// build solution from DG coefficients and the basis
    /// pre-evaluated at a certain point
    Eigen::VectorXd build_sol(const Eigen::VectorXd& u,
                              const Eigen::VectorXd& basis) const;
    
    /// build solution from DG coefficients at a given reference point
    Eigen::VectorXd build_sol(const Eigen::VectorXd& u,
                              double xi) const;

    /// build cell average
    Eigen::MatrixXd build_cell_avg(const Eigen::MatrixXd& u) const;

    /// return boolean if DG limiting is needed
    bool bool_limit_sol() const {
        if (n_coeff == 1) // no limiting for piecewise constant
            return false;
        return true;
    }

    /// build split solution uSol_m = u0 + um, uSol_p = u0 - up
    /// from DG coefficients
    std::tuple <Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
    build_split_sol(const Eigen::MatrixXd& u) const;

    /// build DG coefficients from uSol_m = u0 + um, uSol_p = u0 - up
    void compute_limit_coeffs (Eigen::MatrixXd &u,
                               Eigen::MatrixXd &um,
                               Eigen::MatrixXd &up) const;

  private:
    int n_vars, n_coeff;
    PolynomialBasis poly_basis;
    
    mutable int n_quad;
    mutable Eigen::VectorXd quad_points, quad_weights;
};


#endif // HYPSYS1D_DG_HANDLER_HPP
