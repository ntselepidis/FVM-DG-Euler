#include <ancse/model.hpp>

#include <iostream>
#include <fmt/format.h>


///------------------///
/// Euler equations  ///
///------------------///
Eigen::VectorXd Euler::flux(const Eigen::VectorXd &u) const
{
    return Eigen::VectorXd::Zero(n_vars);
}

Eigen::VectorXd Euler::eigenvalues(const Eigen::VectorXd &u) const
{
    return Eigen::VectorXd::Zero(n_vars);
}

Eigen::MatrixXd Euler::eigenvectors(const Eigen::VectorXd &u) const
{
    return Eigen::MatrixXd::Zero(n_vars, n_vars);
}

double Euler::max_eigenvalue(const Eigen::VectorXd &u) const
{
    return 0;
}


Eigen::VectorXd Euler::cons_to_prim(const Eigen::VectorXd &u_cons) const
{
    return Eigen::VectorXd::Zero(n_vars);
}

Eigen::VectorXd Euler::prim_to_cons(const Eigen::VectorXd &u_prim) const
{
    return Eigen::VectorXd::Zero(n_vars);
}


#define REGISTER_MODEL(token, ModelType)      \
    if (config["model"] == (token)) {         \
        return std::make_shared<ModelType>(); \
    }

std::shared_ptr<Model> make_model (const nlohmann::json &config)
{
    REGISTER_MODEL("burgers", Burgers)

    // implement and register your models here

    throw std::runtime_error(
        fmt::format("Unknown model. {}", std::string(config["flux"])));
}
