#ifndef LBM_POISEUILLE_FLOW_SIMULATOR_HPP
#define LBM_POISEUILLE_FLOW_SIMULATOR_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <chrono>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

struct Parameters {
  std::array<int, 2> grid_shape;
  std::array<double, 2> external_force;
  double relaxation_time;
  double error_limit;
  int print_frequency;
  int max_iter;
};

class PoiseuilleFlowSimulator {
 public:
  PoiseuilleFlowSimulator(const Parameters& params)
      : grid_{params.grid_shape},
        c_{},
        w_{},
        tau_{params.relaxation_time},
        g_{},
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
    g_ = w_.cwiseProduct(c_.transpose() * (3.0 * g));
  }

  Eigen::VectorXd calc_velocity() const noexcept {
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
      tsteps += 1;
      this->calc_equilibrium_distribution_function(feq, u, rho);
      this->collision_process(f, feq, rho);
      this->propagation_process(f, fold);
      this->apply_boundary_condition(f);

      // Compute properties
      rho = f.colwise().sum().transpose();
      u = c_ * f;
      u.array() /= rho.transpose().replicate<2, 1>().array();

      eps = (u - u0).colwise().norm().maxCoeff();
      u0 = u;

      if (tsteps % print_freq_ == 0) {
        fmt::print("iter = {}, eps = {:.6e}\n", tsteps, eps);
      }
    } while (eps > error_limit_ && tsteps < max_iter_);

    const auto end = system_clock::now();

    // in seconds
    const auto elapsed_time =
        duration_cast<milliseconds>(end - start).count() * 1e-3;

    fmt::print("total iter = {}, eps = {:.6e}, simulation time = {:.3} sec\n",
               tsteps, eps, elapsed_time);

    return this->strip_ux_along_y_axis(u);
  }

  template <typename T1, typename T2, typename T3>
  void calc_equilibrium_distribution_function(
      Eigen::MatrixBase<T1>& feq, const Eigen::MatrixBase<T2>& u,
      const Eigen::MatrixBase<T3>& rho) const noexcept {
    for (int i = 0; i < feq.cols(); ++i) {
      const auto u2 = u.col(i).squaredNorm();
      for (int k = 0; k < 9; ++k) {
        const auto cu = c_.col(k).dot(u.col(i));
        feq(k, i) =
            w_(k) * rho(i) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void collision_process(Eigen::MatrixBase<T1>& f,
                         const Eigen::MatrixBase<T2>& feq,
                         const Eigen::MatrixBase<T3>& rho) const noexcept {
    f = f - (f - feq) / tau_ + g_ * rho.transpose();
  }

  template <typename T1, typename T2>
  void propagation_process(Eigen::MatrixBase<T1>& f,
                           Eigen::MatrixBase<T2>& fold) const noexcept {
    fold = f;
    for (int i = 0; i < grid_.ni(); ++i) {
      for (int j = 0; j < grid_.nj(); ++j) {
        const auto n = grid_.index(i, j);
        f(0, n) = fold(0, n);
        f(1, grid_.periodic_index(i, j + 1)) = fold(1, n);
        f(2, grid_.periodic_index(i + 1, j)) = fold(2, n);
        f(3, grid_.periodic_index(i, j - 1)) = fold(3, n);
        f(4, grid_.periodic_index(i - 1, j)) = fold(4, n);
        f(5, grid_.periodic_index(i + 1, j + 1)) = fold(5, n);
        f(6, grid_.periodic_index(i + 1, j - 1)) = fold(6, n);
        f(7, grid_.periodic_index(i - 1, j - 1)) = fold(7, n);
        f(8, grid_.periodic_index(i - 1, j + 1)) = fold(8, n);
      }
    }
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    // Bounce-back condition for bottom boundary
    for (int j = 0; j < grid_.nj(); ++j) {
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
  }

  template <typename T>
  Eigen::VectorXd strip_ux_along_y_axis(
      const Eigen::MatrixBase<T>& u) const noexcept {
    using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Unaligned,
        Eigen::Dynamic;
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), grid_.nj(), grid_.ni(),
        Stride<Dynamic, 2>(grid_.nj() * 2, 2));
    return ux.transpose().col(grid_.nj() / 2);
  }

 private:
  CartesianGrid2d grid_;
  Eigen::Matrix<double, 2, 9> c_;
  Eigen::Matrix<double, 9, 1> w_;
  Eigen::Matrix<double, 9, 1> g_;
  double tau_;
  double error_limit_;
  int print_freq_;
  int max_iter_;
};

}  // namespace lbm

#endif  // LBM_POISEUILLE_FLOW_SIMULATOR_HPP