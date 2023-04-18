#include "lbm/cavity_flow_simulator.hpp"

#include <chrono>

namespace {

double calc_dynamic_viscosity(double relaxation_time) noexcept {
  return (relaxation_time - 0.5) / 3;
}

double calc_reynolds_number(double dynamic_viscosity, double length,
                            double velocity) noexcept {
  return velocity * length / dynamic_viscosity;
}

}  // anonymous namespace

namespace lbm {

CavityFlowSimulator::CavityFlowSimulator(const CavityFlowParameters& params)
    : grid_{params.grid_shape},
      c_{Lattice<LatticeType::D2Q9>::get_lattice_vector()},
      w_{Lattice<LatticeType::D2Q9>::get_weight()},
      dynamic_viscosity_{
          calc_dynamic_viscosity(get_relaxation_time(params.collision_params))},
      length_{static_cast<double>(params.grid_shape[0] - 2)},
      wall_velocity_{params.wall_velocity},
      reynolds_number_{
          calc_reynolds_number(dynamic_viscosity_, length_, wall_velocity_)},
      error_limit_{params.error_limit},
      print_freq_{params.print_frequency},
      max_iter_{params.max_iter},
      writer_{params.output_directory},
      collision_{create_collision_model(params.collision_params)},
      propagator_{grid_},
      north_{grid_, {params.wall_velocity, 0}},
      south_{grid_, {0, 0}},
      east_{grid_, {0, 0}},
      west_{grid_, {0, 0}} {}

void CavityFlowSimulator::run() const {
  using Eigen::MatrixXd, Eigen::VectorXd;
  using Matrix2Xd = Eigen::Matrix<double, 2, Eigen::Dynamic>;
  using Matrix9Xd = Eigen::Matrix<double, 9, Eigen::Dynamic>;
  using std::chrono::system_clock, std::chrono::duration_cast,
      std::chrono::milliseconds;

  // Initialization
  const auto size = grid_.size();
  Matrix2Xd u = Matrix2Xd::Zero(2, size);
  MatrixXd u0 = u;
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
    this->run_propagation_process(f, fold);
    this->apply_boundary_condition(f);
    calc_density(rho, f);
    this->calc_velocity(u, rho, f);

    eps = (u - u0).colwise().norm().maxCoeff();
    u0 = u;

    if (print_freq_ > 0 && tsteps % print_freq_ == 0) {
      fmt::print("iter = {}, eps = {:.6e}\n", tsteps, eps);
    }
  } while (eps > error_limit_ && tsteps < max_iter_);

  const auto end = system_clock::now();

  // in seconds
  const auto elapsed_time =
      duration_cast<milliseconds>(end - start).count() * 1e-3;

  fmt::print(
      "total iter = {}, eps = {:.6e}, Re = {:.3e}, dynamic viscosity = "
      "{:.3}, elapsed time = {:.3f} sec\n",
      tsteps, eps, reynolds_number_, dynamic_viscosity_, elapsed_time);

  this->write_velocity(u);
  this->write_xy();
}

}  // namespace lbm
