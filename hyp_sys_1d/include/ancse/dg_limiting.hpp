#ifndef HYPSYS1D_DG_LIMITING_HPP
#define HYPSYS1D_DG_LIMITING_HPP

#include <ancse/grid.hpp>
#include <ancse/limiters.hpp>
#include <ancse/dg_handler.hpp>
#include <ancse/boundary_condition.hpp>

/// Interface for DGLimiting.
class Limiting {
  public:
    virtual ~Limiting() = default;

    virtual void operator() (Eigen::MatrixXd &u) const = 0;
};


template <class Limiter>
class DGLimiting : public Limiting {
  public:
    DGLimiting (const Grid &grid,
                const DGHandler &dg_handler,
                const Limiter &limiter)
        : n_cells (grid.n_cells),
          n_ghost (grid.n_ghost),
          dx (grid.dx),
          dg_handler (dg_handler),
          limiter (limiter) {}
    
    /// applies the DG limiting procedure to the solution coefficients
    /// for all cells.
    void operator() (Eigen::MatrixXd &u) const;
    
    /// DG limiting procedure
    /** uc0 : cell average of cell 'i'
     *  um0 : cell average of cell 'i-1'
     *  up0 : cell average of cell 'i+1'
     *  uSol_{@right_face_of_cell_i} = uc0 + um
     *  uSol_{@left_face_of_cell_i}  = uc0 - up
    */
    inline std::pair<Eigen::VectorXd, Eigen::VectorXd>
    operator()(const Eigen::VectorXd &uc0,
               const Eigen::VectorXd &um0,
               const Eigen::VectorXd &up0,
               const Eigen::VectorXd &um,
               const Eigen::VectorXd &up) const;

  private:
    int n_cells, n_ghost;
    double dx;

    DGHandler dg_handler;
    Limiter limiter;
};

std::shared_ptr<Limiting>
make_dg_limiting(const nlohmann::json &config,
                 const Grid &grid,
                 const DGHandler &dg_handler);


#endif // HYPSYS1D_DG_LIMITING_HPP
