#ifndef LBM_HALFWAY_BOUNCE_BACK_BOUNDARY_HPP
#define LBM_HALFWAY_BOUNCE_BACK_BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/boundary_type.hpp"
#include "lbm/cartesian_grid_2d.hpp"
#include "lbm/lattice.hpp"

namespace lbm {

template <BoundaryType B>
class HalfwayBounceBackBoundary {
 public:
  HalfwayBounceBackBoundary(const CartesianGrid2d& grid,
                            const std::array<double, 2>& u)
      : grid_{grid}, u_{u[0], u[1]} {}

  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    using Eigen::Matrix;
    const Matrix<double, 2, 9> c =
        Lattice<LatticeType::D2Q9>::get_lattice_vector();
    const Matrix<double, 9, 1> w = Lattice<LatticeType::D2Q9>::get_weight();
    if constexpr (B == BoundaryType::North) {
      const auto north = grid_.ny() - 1;
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, north - 1);
        const auto m1 = grid_.index(i, north);
        const auto m2 = grid_.index(i + 1, north);
        const auto m3 = grid_.index(i - 1, north);
        const auto rho = f(0, m1) + f(1, m1) + f(3, m1) +
                         2.0 * (f(2, m1) + f(5, m1) + f(6, m1));
        f(4, n) = f(2, m1) - 6 * rho * w(2) * (c.col(2).dot(u_));
        f(7, n) = f(5, m2) - 6 * rho * w(5) * (c.col(5).dot(u_));
        f(8, n) = f(6, m3) - 6 * rho * w(6) * (c.col(6).dot(u_));
      }
    } else if constexpr (B == BoundaryType::South) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, 1);
        const auto m1 = grid_.index(i, 0);
        const auto m2 = grid_.index(i - 1, 0);
        const auto m3 = grid_.index(i + 1, 0);
        const auto rho = f(0, m1) + f(1, m1) + f(3, m1) +
                         2.0 * (f(4, m1) + f(7, m1) + f(8, m1));
        f(2, n) = f(4, m1) - 6 * rho * w(4) * (c.col(4).dot(u_));
        f(5, n) = f(7, m2) - 6 * rho * w(7) * (c.col(7).dot(u_));
        f(6, n) = f(8, m3) - 6 * rho * w(8) * (c.col(8).dot(u_));
      }
    } else if constexpr (B == BoundaryType::East) {
      const auto east = grid_.nx() - 1;
      for (int j = 1; j < grid_.ny() - 1; ++j) {
        const auto n = grid_.index(east - 1, j);
        const auto m1 = grid_.index(east, j);
        const auto m2 = grid_.index(east, j - 1);
        const auto m3 = grid_.index(east, j + 1);
        const auto rho = f(0, m1) + f(2, m1) + f(4, m1) +
                         2.0 * (f(1, m1) + f(5, m1) + f(8, m1));
        f(3, n) = f(1, m1) - 6 * rho * w(1) * (c.col(1).dot(u_));
        f(6, n) = f(8, m2) - 6 * rho * w(8) * (c.col(8).dot(u_));
        f(7, n) = f(5, m3) - 6 * rho * w(5) * (c.col(5).dot(u_));
      }
    } else if constexpr (B == BoundaryType::West) {
      for (int j = 1; j < grid_.ny() - 1; ++j) {
        const auto n = grid_.index(1, j);
        const auto m1 = grid_.index(0, j);
        const auto m2 = grid_.index(0, j + 1);
        const auto m3 = grid_.index(0, j - 1);
        const auto rho = f(0, m1) + f(2, m1) + f(4, m1) +
                         2.0 * (f(3, m1) + f(6, m1) + f(7, m1));
        f(1, n) = f(3, m1) - 6 * rho * w(3) * (c.col(3).dot(u_));
        f(8, n) = f(6, m2) - 6 * rho * w(6) * (c.col(6).dot(u_));
        f(5, n) = f(7, m3) - 6 * rho * w(7) * (c.col(7).dot(u_));
      }
    }
  }

 private:
  CartesianGrid2d grid_;
  Eigen::Vector2d u_;
};

}  // namespace lbm

#endif  // LBM_HALFWAY_BOUNCE_BACK_BOUNDARY_HPP