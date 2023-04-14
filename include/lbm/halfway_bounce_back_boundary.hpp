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
      : grid_{grid}, corr_{} {
    using Eigen::all, Eigen::Map, Eigen::Vector3d, Eigen::Matrix,
        Eigen::Vector2d;
    const auto c = Lattice<LatticeType::D2Q9>::get_lattice_vector();
    const auto w = Lattice<LatticeType::D2Q9>::get_weight();
    Matrix<double, 2, 3> ct;
    Vector3d wt;
    if constexpr (B == BoundaryType::North) {
      ct = c(all, {2, 5, 6});
      wt = w({2, 5, 6});
    } else if constexpr (B == BoundaryType::South) {
      ct = c(all, {4, 7, 8});
      wt = w({4, 7, 8});
    } else if constexpr (B == BoundaryType::East) {
      ct = c(all, {1, 8, 5});
      wt = w({1, 8, 5});
    } else if constexpr (B == BoundaryType::West) {
      ct = c(all, {3, 6, 7});
      wt = w({3, 6, 7});
    }
    Map<const Vector2d> ut(u.data());
    corr_ = (6 * wt.array() * (ct.transpose() * ut).array()).matrix();
  }

  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    if constexpr (B == BoundaryType::North) {
      const auto north = grid_.ny() - 1;
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, north - 1);
        const auto m1 = grid_.index(i, north);
        const auto m2 = grid_.index(i + 1, north);
        const auto m3 = grid_.index(i - 1, north);
        const auto rho = f({0, 1, 3}, m1).sum() + 2.0 * f({2, 5, 6}, m1).sum();
        f(4, n) = f(2, m1) - rho * corr_(0);
        f(7, n) = f(5, m2) - rho * corr_(1);
        f(8, n) = f(6, m3) - rho * corr_(2);
      }
    } else if constexpr (B == BoundaryType::South) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, 1);
        const auto m1 = grid_.index(i, 0);
        const auto m2 = grid_.index(i - 1, 0);
        const auto m3 = grid_.index(i + 1, 0);
        const auto rho = f({0, 1, 3}, m1).sum() + 2.0 * f({4, 7, 8}, m1).sum();
        f(2, n) = f(4, m1) - rho * corr_(0);
        f(5, n) = f(7, m2) - rho * corr_(1);
        f(6, n) = f(8, m3) - rho * corr_(2);
      }
    } else if constexpr (B == BoundaryType::East) {
      const auto east = grid_.nx() - 1;
      for (int j = 1; j < grid_.ny() - 1; ++j) {
        const auto n = grid_.index(east - 1, j);
        const auto m1 = grid_.index(east, j);
        const auto m2 = grid_.index(east, j - 1);
        const auto m3 = grid_.index(east, j + 1);
        const auto rho = f({0, 2, 4}, m1).sum() + 2.0 * f({1, 5, 8}, m1).sum();
        f(3, n) = f(1, m1) - rho * corr_(0);
        f(6, n) = f(8, m2) - rho * corr_(1);
        f(7, n) = f(5, m3) - rho * corr_(2);
      }
    } else if constexpr (B == BoundaryType::West) {
      for (int j = 1; j < grid_.ny() - 1; ++j) {
        const auto n = grid_.index(1, j);
        const auto m1 = grid_.index(0, j);
        const auto m2 = grid_.index(0, j + 1);
        const auto m3 = grid_.index(0, j - 1);
        const auto rho = f({0, 2, 4}, m1).sum() + 2.0 * f({3, 6, 7}, m1).sum();
        f(1, n) = f(3, m1) - rho * corr_(0);
        f(8, n) = f(6, m2) - rho * corr_(1);
        f(5, n) = f(7, m3) - rho * corr_(2);
      }
    }
  }

 private:
  CartesianGrid2d grid_;
  Eigen::Vector3d corr_;
};

}  // namespace lbm

#endif  // LBM_HALFWAY_BOUNCE_BACK_BOUNDARY_HPP