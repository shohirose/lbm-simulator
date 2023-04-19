#ifndef LBM_CAVITY_FLOW_SIMULATOR_HPP
#define LBM_CAVITY_FLOW_SIMULATOR_HPP

#include <array>
#include <cassert>
#include <string>

#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/collision_model.hpp"
#include "lbm/file_writer.hpp"
#include "lbm/halfway_bounce_back_boundary.hpp"
#include "lbm/lattice.hpp"
#include "lbm/propagator.hpp"

namespace lbm {

struct CavityFlowParameters {
  std::array<int, 2> grid_shape;  ///< Grid shape [nx, ny]
  double wall_velocity;
  double error_limit;
  int print_frequency;
  int max_iter;
  std::string output_directory;
  CollisionParameters collision_params;
};

class CavityFlowSimulator {
 public:
  using Parameters = CavityFlowParameters;

  /**
   * @brief Construct a new CavityFlowSimulator object
   *
   * @param params Parameters
   * @throw std::filesystem::filesystem_error
   */
  CavityFlowSimulator(const CavityFlowParameters& params);

  /**
   * @brief Run a simulation.
   *
   * @throw std::runtime_error
   */
  void run() const;

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

  template <typename T1, typename T2, typename T3>
  void run_collision_process(Eigen::MatrixBase<T1>& f,
                             const Eigen::MatrixBase<T2>& feq,
                             const Eigen::MatrixBase<T3>& u) const {
    if (std::holds_alternative<SingleRelaxationTimeModel>(collision_)) {
      std::get<SingleRelaxationTimeModel>(collision_).apply(f, feq);
    } else if (std::holds_alternative<MultipleRelaxationTimeModel>(
                   collision_)) {
      std::get<MultipleRelaxationTimeModel>(collision_).apply(f, feq);
    } else if (std::holds_alternative<CentralMomentModel>(
                   collision_)) {
      std::get<CentralMomentModel>(collision_).apply(f, feq, u);
    } else {
      throw std::runtime_error(
          "Error: CollisionModel holds neither of "
          "SingleRelaxationTimeModel, MultipleRelaxationTimeModel, or "
          "CentralMomentModel.");
    }
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
  double dynamic_viscosity_;
  double length_;
  double wall_velocity_;
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