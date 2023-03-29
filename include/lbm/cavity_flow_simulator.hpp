#ifndef LBM_CAVITY_FLOW_SIMULATOR_HPP
#define LBM_CAVITY_FLOW_SIMULATOR_HPP

#include <array>
#include <cassert>
#include <chrono>
#include <nlohmann/json.hpp>
#include <string>

#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/file_writer.hpp"

namespace lbm {

struct CavityFlowParameters {
  std::array<int, 2> grid_shape;
  double wall_velocity;
  double reynolds_number;
  double error_limit;
  int print_frequency;
  int max_iter;
  std::string output_directory;
};

void from_json(const nlohmann::json& j, CavityFlowParameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("wallVelocity").get_to(params.wall_velocity);
  j.at("reynoldsNumber").get_to(params.reynolds_number);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
  j.at("outputDirectory").get_to(params.output_directory);
}

class CavityFlowSimulator {
 public:
  CavityFlowSimulator(const CavityFlowParameters& params)
      : grid_{params.grid_shape},
        c_{},
        w_{},
        ux_{params.wall_velocity},
        tau_{calc_tau(params.reynolds_number, params.wall_velocity,
                      static_cast<double>(grid_.nj() - 1))},
        error_limit_{params.error_limit},
        print_freq_{params.print_frequency},
        max_iter_{params.max_iter},
        writer_{params.output_directory} {
    this->setup();
  }

  CavityFlowSimulator(const std::filesystem::path& p)
      : grid_{},
        c_{},
        w_{},
        ux_{},
        tau_{},
        error_limit_{},
        print_freq_{},
        max_iter_{},
        writer_{} {
    std::ifstream file(p.string());
    if (!file) {
      fmt::print(stderr, "Error: could not open a file: {}", p.string());
      std::exit(EXIT_FAILURE);
    }
    try {
      const auto j = nlohmann::json::parse(file);
      const auto params = j.get<CavityFlowParameters>();
      grid_ = params.grid_shape;
      ux_ = params.wall_velocity;
      tau_ = calc_tau(params.reynolds_number, params.wall_velocity,
                      static_cast<double>(grid_.nj() - 1));
      error_limit_ = params.error_limit;
      print_freq_ = params.print_frequency;
      max_iter_ = params.max_iter;
      writer_.set_output_directory(params.output_directory);
      this->setup();
    } catch (const std::exception& e) {
      fmt::print(stderr, "Error: could not read a JSON file: {}\n", p.string());
      throw e;
    }
  }

  void run() const noexcept {
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
      this->collision_process(f, feq);
      this->propagation_process(f, fold);
      this->apply_boundary_condition(f);
      this->calc_properties(u, rho, f);

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
  static double calc_tau(double Re, double U, double L) {
    const auto nu = U * L / Re;
    return 3.0 * nu + 0.5;
  }

  void setup() {
    // clang-format off
    c_ << 0,  1,  0, -1,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1;
    w_ << 4.0 / 9.0,
          1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
    // clang-format on
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

  template <typename T1, typename T2>
  void collision_process(Eigen::MatrixBase<T1>& f,
                         const Eigen::MatrixBase<T2>& feq) const noexcept {
    f = f - (f - feq) / tau_;
  }

  template <typename T1, typename T2>
  void propagation_process(Eigen::MatrixBase<T1>& f,
                           Eigen::MatrixBase<T2>& fold) const noexcept {
    fold = f;
    for (int i = 1; i < grid_.ni() - 1; ++i) {
      for (int j = 1; j < grid_.nj() - 1; ++j) {
        const auto n = grid_.index(i, j);
        f(1, grid_.index(i, j + 1)) = fold(1, n);
        f(2, grid_.index(i + 1, j)) = fold(2, n);
        f(3, grid_.index(i, j - 1)) = fold(3, n);
        f(4, grid_.index(i - 1, j)) = fold(4, n);
        f(5, grid_.index(i + 1, j + 1)) = fold(5, n);
        f(6, grid_.index(i + 1, j - 1)) = fold(6, n);
        f(7, grid_.index(i - 1, j - 1)) = fold(7, n);
        f(8, grid_.index(i - 1, j + 1)) = fold(8, n);
      }
    }
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    // Half-way bounce-back condition
    const auto right = grid_.nj() - 1;
    const auto top = grid_.ni() - 1;
    // Left boundary
    for (int i = 1; i < grid_.ni() - 1; ++i) {
      const auto n = grid_.index(i, 1);
      f(1, n) = f(3, grid_.index(i, 0));
      f(5, n) = f(7, grid_.index(i - 1, 0));
      f(8, n) = f(6, grid_.index(i + 1, 0));
    }
    // Right boundary
    for (int i = 1; i < grid_.ni() - 1; ++i) {
      const auto n = grid_.index(i, right - 1);
      f(3, n) = f(1, grid_.index(i, right));
      f(6, n) = f(8, grid_.index(i + 1, right));
      f(7, n) = f(5, grid_.index(i - 1, right));
    }
    // Bottom boundary
    for (int j = 1; j < grid_.nj() - 1; ++j) {
      const auto n = grid_.index(1, j);
      f(2, n) = f(4, grid_.index(0, j));
      f(5, n) = f(7, grid_.index(0, j - 1));
      f(6, n) = f(8, grid_.index(0, j + 1));
    }
    // Top boundary
    for (int j = 1; j < grid_.nj() - 1; ++j) {
      const auto n = grid_.index(top - 1, j);
      const auto m = grid_.index(top, j);
      const auto rho =
          f(0, m) + f(1, m) + f(3, m) + 2.0 * (f(2, m) + f(5, m) + f(6, m));
      f(4, n) = f(2, grid_.index(top, j));
      f(7, n) = f(5, grid_.index(top, j + 1)) - rho * ux_ / 6.0;
      f(8, n) = f(6, grid_.index(top, j - 1)) + rho * ux_ / 6.0;
    }
    // Bounce back from corner points
    // f is calculated by using the equation for feq with density = 1.
    const auto u2 = ux_ * ux_;
    {
      const auto cu = c_.col(5).dot(Eigen::Vector2d(ux_, 0.0));
      f(5, grid_.index(top - 1, right - 1)) =
          (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2) / 36.0;
    }
    {
      const auto cu = c_.col(7).dot(Eigen::Vector2d(ux_, 0.0));
      f(7, grid_.index(top - 1, right - 1)) =
          (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2) / 36.0;
    }
    {
      const auto cu = c_.col(6).dot(Eigen::Vector2d(ux_, 0.0));
      f(6, grid_.index(top - 1, 1)) =
          (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2) / 36.0;
    }
    {
      const auto cu = c_.col(8).dot(Eigen::Vector2d(ux_, 0.0));
      f(8, grid_.index(top - 1, 1)) =
          (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2) / 36.0;
    }
    f(6, grid_.index(1, top - 1)) = 1.0 / 36.0;
    f(8, grid_.index(1, top - 1)) = 1.0 / 36.0;
    f(5, grid_.index(1, 1)) = 1.0 / 36.0;
    f(7, grid_.index(1, 1)) = 1.0 / 36.0;
  }

  template <typename T1, typename T2, typename T3>
  void calc_properties(Eigen::MatrixBase<T1>& u, Eigen::MatrixBase<T2>& rho,
                       const Eigen::MatrixBase<T3>& f) const noexcept {
    rho = f.colwise().sum().transpose();
    for (int i = 0; i < u.cols(); ++i) {
      u.col(i) = (c_ * f.col(i)) / rho(i);
    }
  }

  template <typename T>
  void write_velocity(const Eigen::MatrixBase<T>& u) const noexcept {
    using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Unaligned,
        Eigen::Dynamic;
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), grid_.nj(), grid_.ni(),
        Stride<Dynamic, 2>(grid_.nj() * 2, 2));
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> uy(
        &u(1, 0), grid_.nj(), grid_.ni(),
        Stride<Dynamic, 2>(grid_.nj() * 2, 2));
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