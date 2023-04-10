#ifndef LBM_PERIODIC_BOUNDARY_HPP
#define LBM_PERIODIC_BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

enum class PeriodicType {
  NorthSouth,
  EastWest,
  TopBottom,
};

template <PeriodicType P>
struct PeriodicIndexPair {};

template <>
struct PeriodicIndexPair<PeriodicType::NorthSouth> {
  int north;
  int south;

  PeriodicIndexPair(int north_, int south_) : north{north_}, south{south_} {}
};

template <>
struct PeriodicIndexPair<PeriodicType::EastWest> {
  int east;
  int west;

  PeriodicIndexPair(int east_, int west_) : east{east_}, west{west_} {}
};

template <PeriodicType P>
class PeriodicBoundary {
 public:
  using IndexPair = PeriodicIndexPair<P>;

  /**
   * @brief Construct a new PeriodicBoundary object
   *
   * @param cells Cell index pairs of periodic boundaries
   */
  PeriodicBoundary(const std::vector<IndexPair>& cells) : cells_{cells} {}

  /**
   * @brief Construct a new PeriodicBoundary object
   *
   * @param grid Grid
   */
  PeriodicBoundary(const CartesianGrid2d& grid) : cells_{} {
    if constexpr (P == PeriodicType::NorthSouth) {
      cells_.reserve(grid.nx());
      for (int i = 0; i < grid.nx(); ++i) {
        cells_.emplace_back(grid.index(i, grid.ny() - 1), grid.index(i, 0));
      }
    } else if constexpr (P == PeriodicType::EastWest) {
      cells_.reserve(grid.ny());
      for (int j = 0; j < grid.ny(); ++j) {
        cells_.emplace_back(grid.index(grid.nx() - 1, j), grid.index(0, j));
      }
    }
  }

  /**
   * @brief Apply boundary condition
   *
   * @tparam T
   * @param f Distribution function
   */
  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    if constexpr (P == PeriodicType::NorthSouth) {
      for (auto [north, south] : top_) {
        // clang-format off
        f(2, south) = f(2, north);
        f(4, north) = f(4, south);
        f(5, south) = f(5, north);
        f(6, south) = f(6, north);
        f(7, north) = f(7, south);
        f(8, north) = f(8, south);
        // clang-format on
      }
    } else if constexpr (P == PeriodicType ::EastWest) {
      for (auto [east, west] : cells_) {
        // clang-format off
        f(1, west) = f(1, east);
        f(3, east) = f(3, west);
        f(5, west) = f(5, east);
        f(6, east) = f(6, west);
        f(7, east) = f(7, west);
        f(8, west) = f(8, east);
        // clang-format on
      }
    }
  }

 private:
  /// A list of cell index pairs of top & bottom boundaries.
  std::vector<IndexPair> cells_;
};

}  // namespace lbm

#endif  // LBM_PERIODIC_BOUNDARY_HPP