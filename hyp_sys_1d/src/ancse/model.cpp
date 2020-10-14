#include <ancse/model.hpp>

#include <fmt/format.h>
#include <iostream>

///------------------///
/// Euler equations  ///
///------------------///
Eigen::VectorXd Euler::flux(const Eigen::VectorXd &u) const {
  // assume u is given in conservative variables
  Eigen::VectorXd u_prim = cons_to_prim(u);

  const double rho = u_prim(0);
  const double v = u_prim(1);
  const double p = u_prim(2);
  const double E = u(2);

  Eigen::VectorXd fl(3);
  fl << rho * v, rho * v * v + p, (E + p) * v;

  return fl;
}

Eigen::VectorXd Euler::eigenvalues(const Eigen::VectorXd &u) const {
  // assume u is given in conservative variables
  Eigen::VectorXd u_prim = cons_to_prim(u);

  const double rho = u_prim(0);
  const double v = u_prim(1);
  const double p = u_prim(2);
  const double c = std::sqrt(gamma * p) / rho;

  Eigen::VectorXd eigvals(3);
  eigvals << v - c, v, v + c;

  return eigvals;
}

Eigen::MatrixXd Euler::eigenvectors(const Eigen::VectorXd &u) const {
  // assume u is given in conservative variables
  Eigen::VectorXd u_prim = cons_to_prim(u);

  const double rho = u_prim(0);
  const double v = u_prim(1);
  const double p = u_prim(2);
  const double E = u(2);
  const double H = (E + p) / rho;
  const double c = std::sqrt(gamma * p / rho);

  Eigen::MatrixXd eigvecs(3, 3);

  eigvecs.col(0) << 1.0, v - c, H - v * c;
  eigvecs.col(1) << 1.0, v, 0.5 * v * v;
  eigvecs.col(2) << 1.0, v + c, H + v * c;

  return eigvecs;
}

double Euler::max_eigenvalue(const Eigen::VectorXd &u) const {
  return (eigenvalues(u).cwiseAbs()).maxCoeff();
}

Eigen::VectorXd Euler::cons_to_prim(const Eigen::VectorXd &u_cons) const {
  const double rho = u_cons(0);  // density
  const double rhov = u_cons(1); // momentum
  const double E = u_cons(2);    // energy
  const double v = rhov / rho;
  const double p = (E - 0.5 * rhov * v) * (gamma - 1);

  Eigen::VectorXd u_prim(3);
  u_prim << rho, v, p;

  return u_prim;
}

Eigen::VectorXd Euler::prim_to_cons(const Eigen::VectorXd &u_prim) const {
  const double rho = u_prim(0);                         // density
  const double v = u_prim(1);                           // velocity
  const double p = u_prim(2);                           // pressure
  const double E = p / (gamma - 1) + 0.5 * rho * v * v; // energy

  Eigen::VectorXd u_cons(3);
  u_cons << rho, rho * v, E;

  return u_cons;
}

#define REGISTER_MODEL(token, ModelType)                                       \
  if (config["model"] == (token)) {                                            \
    return std::make_shared<ModelType>();                                      \
  }

std::shared_ptr<Model> make_model(const nlohmann::json &config) {
  REGISTER_MODEL("burgers", Burgers)

  REGISTER_MODEL("euler", Euler)

  throw std::runtime_error(
      fmt::format("Unknown model. {}", std::string(config["flux"])));
}
