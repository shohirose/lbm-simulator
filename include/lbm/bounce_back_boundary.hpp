#ifndef LBM_BOUNCE_BACK_BOUNDARY_HPP
#define LBM_BOUNCE_BACK_BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

enum class OuterBoundary {
  North,
  South,
  East,
  West,
  Top,
  Bottom,
};

template <OuterBoundary B>
class BounceBackBoundary {
 public:
  BounceBackBoundary() = default;

  /**
   * @brief Construct a new BounceBackBoundary object
   *
   * @param grid Grid
   */
  BounceBackBoundary(const CartesianGrid2d& grid) : cells_{} {
    if constexpr (B == OuterBoundary::North) {
      cells_.reserve(grid.nx());
      for (int i = 0; i < grid.nx(); ++i) {
        cells_.emplace_back(grid.index(i, grid.ny() - 1));
      }
    } else if constexpr (B == OuterBoundary::South) {
      cells_.reserve(grid.nx());
      for (int i = 0; i < grid.nx(); ++i) {
        cells_.emplace_back(grid.index(i, 0));
      }
    } else if constexpr (B == OuterBoundary::East) {
      cells_.reserve(grid.ny());
      for (int j = 0; j < grid.ny(); ++j) {
        cells_.emplace_back(grid.index(grid.nx() - 1, j));
      }
    } else if constexpr (B == OuterBoundary::West) {
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
    if constexpr (B == OuterBoundary::North) {
      for (auto&& cell : cells_) {
        f(4, cell) = f(2, cell);
        f(7, cell) = f(5, cell);
        f(8, cell) = f(6, cell);
      }
    } else if constexpr (B == OuterBoundary::South) {
      for (auto&& cell : cells_) {
        f(2, cell) = f(4, cell);
        f(5, cell) = f(7, cell);
        f(6, cell) = f(8, cell);
      }
    } else if constexpr (B == OuterBoundary::East) {
      for (auto&& cell : cells_) {
        f(3, cell) = f(1, cell);
        f(6, cell) = f(8, cell);
        f(7, cell) = f(5, cell);
      }
    } else if constexpr (B == OuterBoundary::West) {
      for (auto&& cell : cells_) {
        f(1, cell) = f(3, cell);
        f(8, cell) = f(6, cell);
        f(5, cell) = f(7, cell);
      }
    }
  }

 private:
  std::vector<int> cells_;  ///< A list of cell indices of the boundary.
};

}  // namespace lbm

#endif  // LBM_BOUNCE_BACK_BOUNDARY_HPP