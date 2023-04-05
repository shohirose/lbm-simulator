#ifndef LBM_CAVITY_FLOW_SIMULATOR_HPP
#define LBM_CAVITY_FLOW_SIMULATOR_HPP

#include <array>
#include <cassert>
#include <chrono>
#include <string>

#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/file_writer.hpp"

namespace lbm {

struct CavityFlowParameters {
  std::array<int, 2> grid_shape;  ///< Grid shape [nx, ny]
  double wall_velocity;
  double reynolds_number;
  double error_limit;
  int print_frequency;
  int max_iter;
  std::string output_directory;
};

class CavityFlowSimulator {
 public:
  /**
   * @brief Construct a new CavityFlowSimulator object
   *
   * @param params Parameters
   * @throw std::filesystem::filesystem_error
   */
  CavityFlowSimulator(const CavityFlowParameters& params)
      : grid_{params.grid_shape},
        c_{},
        w_{},
        ux_{params.wall_velocity},
        tau_{calc_relaxation_time(params.reynolds_number, params.wall_velocity,
                                  static_cast<double>(grid_.nx() - 1))},
        error_limit_{params.error_limit},
        print_freq_{params.print_frequency},
        max_iter_{params.max_iter},
        writer_{params.output_directory} {
    // clang-format off
    c_ << 0,  1,  0, -1,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1;
    w_ << 4.0 / 9.0,
          1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
    // clang-format on
  }

  /**
   * @brief Run a simulation.
   *
   * @throw std::runtime_error
   */
  void run() const {
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
      this->run_collision_process(f, feq);
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

    fmt::print("total iter = {}, eps = {:.6e}, elapsed time = {:.3} sec\n",
               tsteps, eps, elapsed_time);

    this->write_velocity(u);
  }

 private:
  static double calc_relaxation_time(double reynolds_number, double velocity,
                                     double length) noexcept {
    // dynamic viscosity
    const auto nu = velocity * length / reynolds_number;
    return 3.0 * nu + 0.5;
  }

  template <typename T1, typename T2, typename T3>
  void calc_equilibrium_distribution_function(
      Eigen::MatrixBase<T1>& feq, const Eigen::MatrixBase<T2>& u,
      const Eigen::MatrixBase<T3>& rho) const noexcept {
    const Eigen::VectorXd u2 = u.colwise().squaredNorm().transpose();
    for (int i = 0; i < feq.cols(); ++i) {
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u.col(i);
      feq.col(i) =
          (rho(i) * w_.array() *
           (1.0 + 3.0 * cu.array() + 4.5 * cu.array().square() - 1.5 * u2(i)))
              .matrix();
    }
  }

  template <typename T1, typename T2>
  void run_collision_process(Eigen::MatrixBase<T1>& f,
                             const Eigen::MatrixBase<T2>& feq) const noexcept {
    f = f - (f - feq) / tau_;
  }

  template <typename T1, typename T2>
  void run_propagation_process(Eigen::MatrixBase<T1>& f,
                               Eigen::MatrixBase<T2>& fold) const noexcept {
    fold = f;
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f(1, grid_.index(i + 1, j    )) = fold(1, n);
        f(2, grid_.index(i,     j + 1)) = fold(2, n);
        f(3, grid_.index(i - 1, j    )) = fold(3, n);
        f(4, grid_.index(i    , j - 1)) = fold(4, n);
        f(5, grid_.index(i + 1, j + 1)) = fold(5, n);
        f(6, grid_.index(i - 1, j + 1)) = fold(6, n);
        f(7, grid_.index(i - 1, j - 1)) = fold(7, n);
        f(8, grid_.index(i + 1, j - 1)) = fold(8, n);
        // clang-format on
      }
    }
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    // Half-way bounce-back condition
    const auto right = grid_.nx() - 1;
    const auto top = grid_.ny() - 1;
    // Left boundary
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      const auto n = grid_.index(1, j);
      f(1, n) = f(3, grid_.index(0, j));
      f(5, n) = f(7, grid_.index(0, j - 1));
      f(8, n) = f(6, grid_.index(0, j + 1));
    }
    // Right boundary
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      const auto n = grid_.index(right - 1, j);
      f(3, n) = f(1, grid_.index(right, j));
      f(6, n) = f(8, grid_.index(right, j - 1));
      f(7, n) = f(5, grid_.index(right, j + 1));
    }
    // Bottom boundary
    for (int i = 1; i < grid_.nx() - 1; ++i) {
      const auto n = grid_.index(i, 1);
      f(2, n) = f(4, grid_.index(i, 0));
      f(5, n) = f(7, grid_.index(i - 1, 0));
      f(6, n) = f(8, grid_.index(i + 1, 0));
    }
    // Top boundary
    for (int i = 1; i < grid_.nx() - 1; ++i) {
      const auto n = grid_.index(i, top - 1);
      const auto m = grid_.index(i, top);
      const auto rho =
          f(0, m) + f(1, m) + f(3, m) + 2.0 * (f(2, m) + f(5, m) + f(6, m));
      f(4, n) = f(2, grid_.index(i, top));
      f(7, n) = f(5, grid_.index(i + 1, top)) - rho * ux_ / 6.0;
      f(8, n) = f(6, grid_.index(i - 1, top)) + rho * ux_ / 6.0;
    }
  }

  template <typename T1, typename T2, typename T3>
  void calc_velocity(Eigen::MatrixBase<T1>& u, const Eigen::MatrixBase<T2>& rho,
                     const Eigen::MatrixBase<T3>& f) const noexcept {
    for (int i = 0; i < u.cols(); ++i) {
      u.col(i) = (c_ * f.col(i)) / rho(i);
    }
  }

  template <typename T1, typename T2>
  static void calc_density(Eigen::MatrixBase<T1>& rho,
                           const Eigen::MatrixBase<T2>& f) noexcept {
    rho = f.colwise().sum().transpose();
  }

  /**
   * @brief Write velocity into files.
   *
   * Velocity components are written into text files "ux.txt" and "uy.txt".
   *
   * @param u Velocity
   * @throw std::runtime_error When failed to open files to write.
   */
  template <typename T>
  void write_velocity(const Eigen::MatrixBase<T>& u) const {
    using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Unaligned,
        Eigen::Dynamic;
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), grid_.nx(), grid_.ny(),
        Stride<Dynamic, 2>(grid_.nx() * 2, 2));
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> uy(
        &u(1, 0), grid_.nx(), grid_.ny(),
        Stride<Dynamic, 2>(grid_.nx() * 2, 2));
    writer_.write(ux.transpose(), "ux.txt");
    writer_.write(uy.transpose(), "uy.txt");
  }

 private:
  CartesianGrid2d grid_;
  Eigen::Matrix<double, 2, 9> c_;
  Eigen::Matrix<double, 9, 1> w_;
  double ux_;
  double tau_;
  double error_limit_;
  int print_freq_;
  int max_iter_;
  FileWriter writer_;
};

}  // namespace lbm

#endif  // LBM_CAVITY_FLOW_SIMULATOR_HPP