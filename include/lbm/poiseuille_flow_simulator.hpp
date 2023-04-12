#ifndef LBM_POISEUILLE_FLOW_SIMULATOR_HPP
#define LBM_POISEUILLE_FLOW_SIMULATOR_HPP

#include <array>
#include <cassert>
#include <chrono>
#include <string>

#include "lbm/bounce_back_boundary.hpp"
#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/file_writer.hpp"
#include "lbm/lattice.hpp"
#include "lbm/periodic_boundary.hpp"
#include "lbm/propagator.hpp"

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

class PoiseuilleFlowSimulator {
 public:
  /**
   * @brief Construct a new PoiseuilleFlowSimulator object
   *
   * @param params Parameters
   * @throw std::filesystem::filesystem_error
   */
  PoiseuilleFlowSimulator(const PoiseuilleFlowParameters& params)
      : grid_{params.grid_shape},
        c_{Lattice<LatticeType::D2Q9>::get_lattice_vector()},
        w_{Lattice<LatticeType::D2Q9>::get_weight()},
        g_{},
        tau_{params.relaxation_time},
        error_limit_{params.error_limit},
        print_freq_{params.print_frequency},
        max_iter_{params.max_iter},
        writer_{params.output_directory},
        propagator_{grid_},
        south_{grid_},
        north_{grid_},
        east_west_{grid_} {
    Eigen::Map<const Eigen::Vector2d> g(params.external_force.data());
    g_ = w_.cwiseProduct(c_.transpose() * (3.0 * g));
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
      ++tsteps;
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
    this->write_y();
  }

 private:
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
    propagator_.apply(f, fold);
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    east_west_.apply(f);
    south_.apply(f);
    north_.apply(f);
  }

  template <typename T1, typename T2, typename T3>
  void calc_velocity(Eigen::MatrixBase<T1>& u, const Eigen::MatrixBase<T2>& rho,
                     const Eigen::MatrixBase<T3>& f) const noexcept {
    u.noalias() = c_ * f;
    u.array() /= rho.transpose().replicate<2, 1>().array();
  }

  template <typename T1, typename T2>
  static void calc_density(Eigen::MatrixBase<T1>& rho,
                           const Eigen::MatrixBase<T2>& f) noexcept {
    rho = f.colwise().sum().transpose();
  }

  /**
   * @brief Write the x component of velocity into "ux.txt".
   *
   * @param u Velocity
   * @throw std::runtime_error When failed to open a file to write.
   */
  template <typename T>
  void write_velocity(const Eigen::MatrixBase<T>& u) const {
    using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Unaligned,
        Eigen::Dynamic;
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), grid_.nx(), grid_.ny(),
        Stride<Dynamic, 2>(grid_.nx() * 2, 2));
    writer_.write(ux.transpose().col(grid_.nx() / 2), "ux.txt");
  }

  /**
   * @brief Write the y coordinates into "y.txt".
   *
   * @throw std::runtime_error When failed to open a file to write.
   */
  void write_y() const {
    Eigen::MatrixXd y =
        Eigen::VectorXd::LinSpaced(grid_.ny(), 0, grid_.ny() - 1);
    writer_.write(y, "y.txt");
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
  AllCellPropagator propagator_;
  BounceBackBoundary<BoundaryType::South> south_;
  BounceBackBoundary<BoundaryType::North> north_;
  PeriodicBoundary<PeriodicType::EastWest> east_west_;
};

}  // namespace lbm

#endif  // LBM_POISEUILLE_FLOW_SIMULATOR_HPP