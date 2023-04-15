#ifndef LBM_CAVITY_FLOW_SIMULATOR_HPP
#define LBM_CAVITY_FLOW_SIMULATOR_HPP

#include <array>
#include <cassert>
#include <chrono>
#include <string>

#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/file_writer.hpp"
#include "lbm/halfway_bounce_back_boundary.hpp"
#include "lbm/lattice.hpp"
#include "lbm/propagator.hpp"
#include "lbm/type_traits.hpp"

namespace lbm {

template <typename CollisionParameters>
struct CavityFlowParameters {
  std::array<int, 2> grid_shape;  ///< Grid shape [nx, ny]
  double wall_velocity;
  double error_limit;
  int print_frequency;
  int max_iter;
  std::string output_directory;
  CollisionParameters collision_params;
};

template <typename CollisionModel>
class CavityFlowSimulator {
 public:
  using CollisionParameters = get_collision_parameters_t<CollisionModel>;
  using Parameters = CavityFlowParameters<CollisionParameters>;

  /**
   * @brief Construct a new CavityFlowSimulator object
   *
   * @param params Parameters
   * @throw std::filesystem::filesystem_error
   */
  CavityFlowSimulator(const Parameters& params)
      : grid_{params.grid_shape},
        c_{Lattice<LatticeType::D2Q9>::get_lattice_vector()},
        w_{Lattice<LatticeType::D2Q9>::get_weight()},
        reynolds_number_{calc_reynolds_number(params.collision_params.tau,
                                              params.grid_shape[0] - 2,
                                              params.wall_velocity)},
        error_limit_{params.error_limit},
        print_freq_{params.print_frequency},
        max_iter_{params.max_iter},
        writer_{params.output_directory},
        collision_{params.collision_params},
        propagator_{grid_},
        north_{grid_, {params.wall_velocity, 0}},
        south_{grid_, {0, 0}},
        east_{grid_, {0, 0}},
        west_{grid_, {0, 0}} {}

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
        "total iter = {}, eps = {:.6e}, Re = {:.3e}, elapsed time = {:.3f} "
        "sec\n",
        tsteps, eps, reynolds_number_, elapsed_time);

    this->write_velocity(u);
    this->write_xy();
  }

 private:
  static double calc_reynolds_number(double relaxation_time, double length,
                                     double velocity) noexcept {
    // dynamic viscosity
    const auto nu = (relaxation_time - 0.5) / 3;
    return velocity * length / nu;
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
    collision_.apply(f, feq);
  }

  template <typename T1, typename T2>
  void run_propagation_process(Eigen::MatrixBase<T1>& f,
                               Eigen::MatrixBase<T2>& fold) const noexcept {
    propagator_.apply(f, fold);
  }

  template <typename T>
  void apply_boundary_condition(Eigen::MatrixBase<T>& f) const noexcept {
    east_.apply(f);
    west_.apply(f);
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
   * @brief Write velocity into files.
   *
   * Velocity components are written into "ux.txt" and "uy.txt".
   *
   * @param u Velocity
   * @throw std::runtime_error When failed to open files to write.
   */
  template <typename T>
  void write_velocity(const Eigen::MatrixBase<T>& u) const {
    using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Unaligned,
        Eigen::Dynamic, Eigen::seqN;
    const auto nx = grid_.nx();
    const auto ny = grid_.ny();
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), nx, ny, Stride<Dynamic, 2>(nx * 2, 2));
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> uy(
        &u(1, 0), nx, ny, Stride<Dynamic, 2>(nx * 2, 2));
    writer_.write(ux(seqN(1, nx - 2), seqN(1, ny - 2)).transpose(), "ux.txt");
    writer_.write(uy(seqN(1, nx - 2), seqN(1, ny - 2)).transpose(), "uy.txt");
  }

  /**
   * @brief Write the x and y coordinates into "x.txt" and "y.txt".
   *
   * @throw std::runtime_error When failed to open files to write.
   */
  void write_xy() const {
    using Eigen::VectorXd, Eigen::MatrixXd;
    const auto nx = grid_.nx();
    const auto ny = grid_.ny();
    const VectorXd x = VectorXd::LinSpaced(nx - 2, 0.5, nx - 2.5);
    const VectorXd y = VectorXd::LinSpaced(ny - 2, 0.5, ny - 2.5);
    writer_.write(x.replicate(1, ny - 2).transpose(), "x.txt");
    writer_.write(y.replicate(1, nx - 2), "y.txt");
  }

 private:
  CartesianGrid2d grid_;
  Eigen::Matrix<double, 2, 9> c_;
  Eigen::Matrix<double, 9, 1> w_;
  double reynolds_number_;
  double error_limit_;
  int print_freq_;
  int max_iter_;
  FileWriter writer_;
  CollisionModel collision_;
  InternalCellPropagator propagator_;
  HalfwayBounceBackBoundary<BoundaryType::North> north_;
  HalfwayBounceBackBoundary<BoundaryType::South> south_;
  HalfwayBounceBackBoundary<BoundaryType::East> east_;
  HalfwayBounceBackBoundary<BoundaryType::West> west_;
};

}  // namespace lbm

#endif  // LBM_CAVITY_FLOW_SIMULATOR_HPP