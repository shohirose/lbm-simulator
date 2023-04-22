#ifndef LBM_POISEUILLE_FLOW_SIMULATOR_HPP
#define LBM_POISEUILLE_FLOW_SIMULATOR_HPP

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <string>

#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/collision_model.hpp"
#include "lbm/file_writer.hpp"
#include "lbm/lattice.hpp"
#include "lbm/ongrid_bounce_back_boundary.hpp"
#include "lbm/periodic_boundary.hpp"
#include "lbm/propagator.hpp"

namespace lbm {

struct PoiseuilleFlowParameters {
  std::array<int, 2> grid_shape;
  std::array<double, 2> external_force;
  double error_limit;
  int print_frequency;
  /// Frequency of updating the relative change of velocity
  int relative_change_frequency;
  int max_iter;
  std::string output_directory;
  CollisionParameters collision_params;
};

class PoiseuilleFlowSimulator {
 public:
  using Parameters = PoiseuilleFlowParameters;
  /**
   * @brief Construct a new PoiseuilleFlowSimulator object
   *
   * @param params Parameters
   * @throw std::filesystem::filesystem_error, std::runtime_error
   */
  PoiseuilleFlowSimulator(const PoiseuilleFlowParameters& params);

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
    std::visit([&f, &feq, &u](const auto& model) { model.apply(f, feq, u); },
               collision_);
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
  double error_limit_;
  int print_freq_;
  int eps_freq_;
  int max_iter_;
  FileWriter writer_;
  CollisionModel collision_;
  AllCellPropagator propagator_;
  OnGridBounceBackBoundary<BoundaryType::South> south_;
  OnGridBounceBackBoundary<BoundaryType::North> north_;
  PeriodicBoundary<PeriodicType::EastWest> east_west_;
};

}  // namespace lbm

#endif  // LBM_POISEUILLE_FLOW_SIMULATOR_HPP