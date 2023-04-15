#ifndef LBM_ONGRID_BOUNCE_BACK_BOUNDARY_HPP
#define LBM_ONGRID_BOUNCE_BACK_BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/boundary_type.hpp"
#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

template <BoundaryType B>
class OnGridBounceBackBoundary {
 public:
  OnGridBounceBackBoundary() = default;

  /**
   * @brief Construct a new OnGridBounceBackBoundary object
   *
   * @param grid Grid
   */
  OnGridBounceBackBoundary(const CartesianGrid2d& grid) : cells_{} {
    if constexpr (B == BoundaryType::North) {
      cells_.reserve(grid.nx());
      for (int i = 0; i < grid.nx(); ++i) {
        cells_.emplace_back(grid.index(i, grid.ny() - 1));
      }
    } else if constexpr (B == BoundaryType::South) {
      cells_.reserve(grid.nx());
      for (int i = 0; i < grid.nx(); ++i) {
        cells_.emplace_back(grid.index(i, 0));
      }
    } else if constexpr (B == BoundaryType::East) {
      cells_.reserve(grid.ny());
      for (int j = 0; j < grid.ny(); ++j) {
        cells_.emplace_back(grid.index(grid.nx() - 1, j));
      }
    } else if constexpr (B == BoundaryType::West) {
      cells_.reserve(grid.ny());
      for (int j = 0; j < grid.ny(); ++j) {
        cells_.emplace_back(grid.index(0, j));
      }
    }
  }

  /**
   * @brief Apply the boundary condition
   *
   * @tparam T
   * @param f Distribution function
   */
  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    if constexpr (B == BoundaryType::North) {
      f({4, 7, 8}, cells_) = f({2, 5, 6}, cells_);
    } else if constexpr (B == BoundaryType::South) {
      f({2, 5, 6}, cells_) = f({4, 7, 8}, cells_);
    } else if constexpr (B == BoundaryType::East) {
      f({3, 6, 7}, cells_) = f({1, 8, 5}, cells_);
    } else if constexpr (B == BoundaryType::West) {
      f({1, 8, 5}, cells_) = f({3, 6, 7}, cells_);
    }
  }

 private:
  std::vector<int> cells_;  ///< A list of cell indices of the boundary.
};

}  // namespace lbm

#endif  // LBM_ONGRID_BOUNCE_BACK_BOUNDARY_HPP