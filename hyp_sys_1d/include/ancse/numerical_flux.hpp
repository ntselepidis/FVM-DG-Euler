#ifndef HYPSYS1D_NUMERICAL_FLUX_HPP
#define HYPSYS1D_NUMERICAL_FLUX_HPP

#include <memory>
#include <iostream>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model.
 * It is also unconditionally a bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access
    //       to the following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const std::shared_ptr<Model> &model)
        : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    std::shared_ptr<Model> model;
};

/// Rusanov flux.
class Rusanov {
  public:
    explicit Rusanov (const std::shared_ptr<Model> &model)
        : model(model) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL, 
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);
        
        double lambda = std::max(model->max_eigenvalue(uL), 
                                 model->max_eigenvalue(uR));

        return 0.5 * ( (fL + fR) - lambda*(uR - uL) );
    }

  private:
    std::shared_ptr<Model> model;
};

/// LaxFriedrichs flux.
class LaxFriedrichs {
  public:
    // Note: This version is a bit tricky. A numerical flux should be
    //       a function of the two trace values at the interface, i.e. what we
    //       call `uL`, `uR`. However, it requires 'dt' and 'dx'. Therefore,
    //       these need to be made available to the flux. This is one of the
    //       reasons why `SimulationTime`.
    explicit LaxFriedrichs(
        const Grid &grid,
        const std::shared_ptr<Model> &model,
        const std::shared_ptr<SimulationTime> &simulation_time) 
    : grid(grid), 
      model(model), 
      simulation_time(simulation_time) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL, 
                               const Eigen::VectorXd &uR) const
    {
        const double dx = grid.dx;
        const double dt = simulation_time->dt;

        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * ( (fL + fR) - (dx/dt)*(uR - uL) );
    }

  private:
    std::shared_ptr<Model> model;
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
};

/// Roe flux.
class Roe {
  public:
    explicit Roe(const std::shared_ptr<Model> &model) 
        : model(model) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL, 
                               const Eigen::VectorXd &uR) const
    {
        Eigen::VectorXd uL_prim = model->cons_to_prim( uL );
        Eigen::VectorXd uR_prim = model->cons_to_prim( uR );

        const double rhoL = uL_prim(0);
        const double rhoR = uR_prim(0);
        const double   vL = uL_prim(1);
        const double   vR = uR_prim(1);
        const double   pL = uL_prim(2);
        const double   pR = uR_prim(2);
        const double   EL = uL(2);
        const double   ER = uR(2);
        const double   HL = ( EL + pL ) / rhoL;
        const double   HR = ( ER + pR ) / rhoR;

        const double wL = std::sqrt( rhoL );
        const double wR = std::sqrt( rhoR );

        const double rho_hat = 0.5 * ( rhoL + rhoR );
        const double   v_hat = ( wL*vL + wR*vR ) / ( wL + wR );
        const double   H_hat = ( wL*HL + wR*HR ) / ( wL + wR );

        const double gamma = 
    ( std::dynamic_pointer_cast<Euler>(model) )->get_gamma();
        
        const double E_hat = ( rho_hat * H_hat ) / gamma 
    + 0.5 * rho_hat * v_hat * v_hat * ( 1.0 - 1.0 / gamma ); 

        Eigen::VectorXd u_hat(3);
        u_hat << rho_hat, rho_hat*v_hat, E_hat;

        Eigen::VectorXd L = model->eigenvalues( u_hat );
        Eigen::MatrixXd R = model->eigenvectors( u_hat );
        Eigen::VectorXd x = R.partialPivLu().solve( uR - uL );

        // Harten's entropy fix.
        const double eps = 1e-10;
        for (int i=0; i<3; i++)
            if ( L(i) > eps ) 
                L(i) = std::abs(L(i));
            else
                L(i) = (std::abs(L(i)*L(i)) + eps*eps) / (2*eps); 

        x = R * ( L.array() * x.array() ).matrix();

        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5*( (fL + fR) - x );
    }

  private:
    std::shared_ptr<Model> model;
};

#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
