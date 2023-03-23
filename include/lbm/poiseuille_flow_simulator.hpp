#ifndef LBM_POISEUILLE_FLOW_SIMULATOR_HPP
#define LBM_POISEUILLE_FLOW_SIMULATOR_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <chrono>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

struct Params {
  std::array<int, 2> grid_shape;
  std::array<double, 2> external_force;
  double relaxation_time;
  double error_limit;
  int print_frequency;
  int max_iter;
};

class PoiseuilleFlowSimulator {
 public:
  PoiseuilleFlowSimulator(const Params& params)
      : grid_{params.grid_shape},
        c_{},
        w_{},
        tau_{params.relaxation_time},
        wcg_{},
        error_limit_{params.error_limit},
        print_freq_{params.print_frequency},
        max_iter_{params.max_iter} {
    // clang-format off
    c_ << 0,  1,  0, -1,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1;
    w_ << 4.0 / 9.0,
          1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
    // clang-format on

    Eigen::Map<const Eigen::Vector2d> g(params.external_force.data());
    wcg_ = w_.cwiseProduct(c_.transpose() * (3.0 * g));
  }

  Eigen::MatrixXd calc_velocity() const noexcept {
    const auto size = grid_.size();
    using Eigen::MatrixXd, Eigen::VectorXd;
    using std::chrono::system_clock, std::chrono::duration_cast,
        std::chrono::milliseconds;

    // Initialization
    MatrixXd u = MatrixXd::Zero(2, size);
    MatrixXd u0 = u;
    VectorXd rho = VectorXd::Ones(size);
    MatrixXd feq = this->calc_equilibrium_distribution_function(u, rho);
    MatrixXd f = feq;

    double eps = 1.0;
    int tsteps = 0;
    const auto start = system_clock::now();

    do {
      tsteps += 1;
      feq = this->calc_equilibrium_distribution_function(u, rho);
      f = this->collision_process(f, feq, rho);
      f = this->propagation_process(f);
      this->apply_boundary_condition(f);
      rho.noalias() = f.colwise().sum().transpose();
      u.noalias() = (c_ * f).cwiseQuotient(rho.replicate<1, 2>().transpose());
      eps = (u - u0).norm();
      u0 = u;

      if (tsteps % print_freq_ == 0) {
        fmt::print("iter = {}, eps = {:.6e}\n", tsteps, eps);
      }
    } while (eps > error_limit_ && tsteps < max_iter_);

    const auto end = system_clock::now();
    const auto elapsed_time = duration_cast<milliseconds>(end - start).count();

    fmt::print("total iter = {}, eps = {:.6e}, simulation time = {:.3} sec\n",
               tsteps, eps, elapsed_time * 1e-3);

    return u;
  }

  template <typename T1, typename T2>
  Eigen::MatrixXd calc_equilibrium_distribution_function(
      const Eigen::MatrixBase<T1>& u,
      const Eigen::MatrixBase<T2>& rho) const noexcept {
    using Eigen::MatrixXd;
    MatrixXd feq(9, u.cols());
    for (int i = 0; i < u.cols(); ++i) {
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u.col(i);
      feq.col(i).noalias() =
          (w_.array() * rho(i) *
           (1.0 + 3.0 * cu.array() + 4.5 * cu.array().square() -
            1.5 * u.col(i).squaredNorm()))
              .matrix();
    }
    return feq;
  }

  template <typename T1, typename T2, typename T3>
  Eigen::MatrixXd collision_process(
      const Eigen::MatrixBase<T1>& f, const Eigen::MatrixBase<T2>& feq,
      const Eigen::MatrixBase<T3>& rho) const noexcept {
    return f - (f - feq) / tau_ + wcg_ * rho.transpose();
  }

  template <typename T>
  Eigen::MatrixXd propagation_process(
      const Eigen::MatrixBase<T>& f) const noexcept {
    using Eigen::MatrixXd;
    MatrixXd fnew = MatrixXd::Zero(f.rows(), f.cols());
    for (int i = 0; i < grid_.ni(); ++i) {
      for (int j = 0; j < grid_.nj(); ++j) {
        const auto n = grid_.index(i, j);
        fnew(0, n) = f(0, n);
        fnew(1, grid_.periodic_index(i, j + 1)) = f(1, n);
        fnew(2, grid_.periodic_index(i + 1, j)) = f(2, n);
        fnew(3, grid_.periodic_index(i, j - 1)) = f(3, n);
        fnew(4, grid_.periodic_index(i - 1, j)) = f(4, n);
        fnew(5, grid_.periodic_index(i + 1, j + 1)) = f(5, n);
        fnew(6, grid_.periodic_index(i + 1, j - 1)) = f(6, n);
        fnew(7, grid_.periodic_index(i - 1, j - 1)) = f(7, n);
        fnew(8, grid_.periodic_index(i - 1, j + 1)) = f(8, n);
      }
    }
    return fnew;
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    // Bounce-back condition for bottom boundary
    for (int j = 0; j < grid_.nj(); ++j) {
      //const auto n = grid_.index(0, j);
      const auto n = grid_.index(bottom, j);
      f(2, n) = f(4, n);
      f(5, n) = f(7, n);
      f(6, n) = f(8, n);
    }
    // Bounce-back condition for top boundary
    for (int j = 0; j < grid_.nj(); ++j) {
      const auto n = grid_.index(top, j);
      f(4, n) = f(2, n);
      f(7, n) = f(5, n);
      f(8, n) = f(6, n);
    }
    // Periodic condition
    for (int i = 0; i < grid_.ni(); ++i) {
      const auto nl = grid_.index(i, left);
      const auto nr = grid_.index(i, right);
      f(1, nl) = f(1, nr);
      f(3, nr) = f(3, nl);
      f(5, nl) = f(5, nr);
      f(6, nr) = f(6, nl);
      f(7, nr) = f(7, nl);
      f(8, nl) = f(8, nr);
    }
  }

 private:
  CartesianGrid2d grid_;
  Eigen::Matrix<double, 2, 9> c_;
  Eigen::Matrix<double, 9, 1> w_;
  Eigen::Matrix<double, 9, 1> wcg_;
  double tau_;
  double error_limit_;
  int print_freq_;
  int max_iter_;
};

}  // namespace lbm

#endif  // LBM_POISEUILLE_FLOW_SIMULATOR_HPP