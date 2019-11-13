#include <ancse/dg_limiting.hpp>

#include <iostream>
#include <ancse/config.hpp>
#include <ancse/boundary_condition.hpp>


/// applies the DG limiting procedure to the solution coefficients
/// for all cells.
template <class Limiter>
void DGLimiting <Limiter> :: operator()(Eigen::MatrixXd &u) const
{
}

/// DG limiting procedure
/** uc0 : cell average of cell 'i'
 *  um0 : cell average of cell 'i-1'
 *  up0 : cell average of cell 'i+1'
 *  uSol_{@right_face_of_cell_i} = uc0 + um
 *  uSol_{@left_face_of_cell_i}  = uc0 - up
*/
template <class Limiter>
inline std::pair<Eigen::VectorXd, Eigen::VectorXd>
DGLimiting <Limiter> :: operator()(const Eigen::VectorXd &um0,
                                   const Eigen::VectorXd &uc0,
                                   const Eigen::VectorXd &up0,
                                   const Eigen::VectorXd &um,
                                   const Eigen::VectorXd &up) const
{
    Eigen::VectorXd uml(uc0.size());
    Eigen::VectorXd upl(uc0.size());


    return {std::move(uml), std::move(upl)};
}


#define REGISTER_DG_LIMITER(token, LimiterType, limiter)      \
    if (config["dg_limiter"] == token) {                      \
        return std::make_shared<DGLimiting<LimiterType>> (    \
            grid, dg_handler, limiter);   \
    }

std::shared_ptr<Limiting>
make_dg_limiting(const nlohmann::json &config,
                 const Grid &grid,
                 const DGHandler &dg_handler)
{
    REGISTER_DG_LIMITER("vanleer", VanLeer, VanLeer{})

    REGISTER_DG_LIMITER("shu", Shu, Shu(grid.dx))

    throw std::runtime_error(fmt::format(
        "Unknown DG limiter. [{}]", std::string(config["dg_limiter"])));
}

#undef REGISTER_DG_LIMITER
