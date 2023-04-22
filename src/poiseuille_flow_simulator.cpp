#include "lbm/poiseuille_flow_simulator.hpp"

#include <chrono>

namespace lbm {

PoiseuilleFlowSimulator::PoiseuilleFlowSimulator(
    const PoiseuilleFlowParameters& params)
    : grid_{params.grid_shape},
      c_{Lattice<LatticeType::D2Q9>::get_lattice_vector()},
      w_{Lattice<LatticeType::D2Q9>::get_weight()},
      g_{},
      error_limit_{params.error_limit},
      print_freq_{params.print_frequency},
      eps_freq_{params.relative_change_frequency},
      max_iter_{params.max_iter},
      writer_{params.output_directory},
      collision_{create_collision_model(params.collision_params)},
      propagator_{grid_},
      south_{grid_},
      north_{grid_},
      east_west_{grid_} {
  Eigen::Map<const Eigen::Vector2d> g(params.external_force.data());
  g_ = w_.cwiseProduct(c_.transpose() * (3.0 * g));
}

void PoiseuilleFlowSimulator::run() const {
  using Eigen::MatrixXd, Eigen::VectorXd;
  using Matrix2Xd = Eigen::Matrix<double, 2, Eigen::Dynamic>;
  using Matrix9Xd = Eigen::Matrix<double, 9, Eigen::Dynamic>;
  using std::chrono::system_clock, std::chrono::duration_cast,
      std::chrono::milliseconds;

  // Initialization
  const auto size = grid_.size();
  Matrix2Xd u = Matrix2Xd::Zero(2, size);
  MatrixXd uold = u;
  VectorXd rho = VectorXd::Ones(size);
  Matrix9Xd feq(9, size);
  this->calc_equilibrium_distribution_function(feq, u, rho);
  Matrix9Xd f = feq;
  Matrix9Xd fold = f;

  double eps = 1.0;
  int tsteps = 0;
  const auto start = system_clock::now();

  do {
    ++tsteps;
    this->calc_equilibrium_distribution_function(feq, u, rho);
    this->run_collision_process(f, feq, u);
    this->add_external_force(f, rho);
    this->run_propagation_process(f, fold);
    this->apply_boundary_condition(f);
    calc_density(rho, f);
    this->calc_velocity(u, rho, f);

    if (tsteps % eps_freq_ == 0) {
      eps = ((u - uold).colwise().norm() /
             u.colwise().norm().maxCoeff<Eigen::PropagateNumbers>())
                .maxCoeff();
      uold = u;
    }

    if (print_freq_ > 0 && tsteps % print_freq_ == 0) {
      fmt::print("iter = {}, eps = {:.6e}\n", tsteps, eps);
    }
  } while (eps > error_limit_ && tsteps < max_iter_);

  const auto end = system_clock::now();

  // in seconds
  const auto elapsed_time =
      duration_cast<milliseconds>(end - start).count() * 1e-3;

  fmt::print("total iter = {}, eps = {:.6e}, elapsed time = {:.3} sec\n",
             tsteps, eps, elapsed_time);

  this->write_velocity(u);
  this->write_y();
}

}  // namespace lbm
