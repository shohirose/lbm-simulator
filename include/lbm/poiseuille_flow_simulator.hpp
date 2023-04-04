#ifndef LBM_POISEUILLE_FLOW_SIMULATOR_HPP
#define LBM_POISEUILLE_FLOW_SIMULATOR_HPP

#include <array>
#include <cassert>
#include <chrono>
#include <nlohmann/json.hpp>
#include <string>

#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/file_writer.hpp"

namespace lbm {

struct PoiseuilleFlowParameters {
  std::array<int, 2> grid_shape;
  std::array<double, 2> external_force;
  double relaxation_time;
  double error_limit;
  int print_frequency;
  int max_iter;
  std::string output_directory;
};

void from_json(const nlohmann::json& j, PoiseuilleFlowParameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("externalForce").get_to(params.external_force);
  j.at("relaxationTime").get_to(params.relaxation_time);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
  j.at("outputDirectory").get_to(params.output_directory);
}

class PoiseuilleFlowSimulator {
 public:
  /**
   * @brief Construct a new Poiseuille Flow Simulator object
   *
   * @param params Parameters
   * @throw std::filesystem::filesystem_error
   */
  PoiseuilleFlowSimulator(const PoiseuilleFlowParameters& params)
      : grid_{params.grid_shape},
        c_{},
        w_{},
        g_{},
        tau_{params.relaxation_time},
        error_limit_{params.error_limit},
        print_freq_{params.print_frequency},
        max_iter_{params.max_iter},
        writer_{params.output_directory} {
    this->setup(params.external_force);
  }

  PoiseuilleFlowSimulator(const std::filesystem::path& p)
      : grid_{},
        c_{},
        w_{},
        g_{},
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
      const auto params = j.get<PoiseuilleFlowParameters>();
      grid_ = params.grid_shape;
      tau_ = params.relaxation_time;
      error_limit_ = params.error_limit;
      print_freq_ = params.print_frequency;
      max_iter_ = params.max_iter;
      writer_.set_output_directory(params.output_directory);
      this->setup(params.external_force);
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
      this->run_collision_process(f, feq);
      this->add_external_force(f, rho);
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
  void setup(const std::array<double, 2>& external_force) {
    // clang-format off
    c_ << 0,  1,  0, -1,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1;
    w_ << 4.0 / 9.0,
          1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
    // clang-format on

    Eigen::Map<const Eigen::Vector2d> g(external_force.data());
    g_ = w_.cwiseProduct(c_.transpose() * (3.0 * g));
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
  void add_external_force(Eigen::MatrixBase<T1>& f,
                          const Eigen::MatrixBase<T2>& rho) const noexcept {
    f += g_ * rho.transpose();
  }

  template <typename T1, typename T2>
  void run_propagation_process(Eigen::MatrixBase<T1>& f,
                               Eigen::MatrixBase<T2>& fold) const noexcept {
    fold = f;
    for (int i = 0; i < grid_.nx(); ++i) {
      const auto n = grid_.index(i, 0);
      // clang-format off
      f(1, grid_.periodic_index(i + 1, 0)) = fold(1, n);
      f(2, grid_.periodic_index(i    , 1)) = fold(2, n);
      f(3, grid_.periodic_index(i - 1, 0)) = fold(3, n);
      f(5, grid_.periodic_index(i + 1, 1)) = fold(5, n);
      f(6, grid_.periodic_index(i - 1, 1)) = fold(6, n);
      // clang-format on
    }
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 0; i < grid_.nx(); ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f(1, grid_.periodic_index(i + 1, j    )) = fold(1, n);
        f(2, grid_.periodic_index(i    , j + 1)) = fold(2, n);
        f(3, grid_.periodic_index(i - 1, j    )) = fold(3, n);
        f(4, grid_.periodic_index(i    , j - 1)) = fold(4, n);
        f(5, grid_.periodic_index(i + 1, j + 1)) = fold(5, n);
        f(6, grid_.periodic_index(i - 1, j + 1)) = fold(6, n);
        f(7, grid_.periodic_index(i - 1, j - 1)) = fold(7, n);
        f(8, grid_.periodic_index(i + 1, j - 1)) = fold(8, n);
        // clang-format on
      }
    }
    const auto top = grid_.ny() - 1;
    for (int i = 0; i < grid_.nx() - 1; ++i) {
      const auto n = grid_.index(i, top);
      // clang-format off
      f(1, grid_.periodic_index(i + 1, top    )) = fold(1, n);
      f(3, grid_.periodic_index(i - 1, top    )) = fold(3, n);
      f(4, grid_.periodic_index(i    , top - 1)) = fold(4, n);
      f(7, grid_.periodic_index(i - 1, top - 1)) = fold(7, n);
      f(8, grid_.periodic_index(i + 1, top - 1)) = fold(8, n);
      // clang-format on
    }
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    // Bounce-back condition for bottom boundary
    for (int i = 0; i < grid_.nx(); ++i) {
      const auto n = grid_.index(i, 0);
      f(2, n) = f(4, n);
      f(5, n) = f(7, n);
      f(6, n) = f(8, n);
    }
    // Bounce-back condition for top boundary
    const auto top = grid_.ny() - 1;
    for (int i = 0; i < grid_.nx(); ++i) {
      const auto n = grid_.index(i, top);
      f(4, n) = f(2, n);
      f(7, n) = f(5, n);
      f(8, n) = f(6, n);
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

  template <typename T>
  void write_velocity(const Eigen::MatrixBase<T>& u) const noexcept {
    using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Unaligned,
        Eigen::Dynamic;
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), grid_.nx(), grid_.ny(),
        Stride<Dynamic, 2>(grid_.nx() * 2, 2));
    writer_.write(ux.transpose().col(grid_.nx() / 2), "ux.txt");
  }

  CartesianGrid2d grid_;
  Eigen::Matrix<double, 2, 9> c_;
  Eigen::Matrix<double, 9, 1> w_;
  Eigen::Matrix<double, 9, 1> g_;
  double tau_;
  double error_limit_;
  int print_freq_;
  int max_iter_;
  FileWriter writer_;
};

}  // namespace lbm

#endif  // LBM_POISEUILLE_FLOW_SIMULATOR_HPP