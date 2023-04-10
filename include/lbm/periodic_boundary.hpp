#ifndef LBM_PERIODIC_BOUNDARY_HPP
#define LBM_PERIODIC_BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

struct TopBottomIndexPair {
  int top;     ///< Cell index of top boundary
  int bottom;  ///< Cell index of bottom boundary

  TopBottomIndexPair(int top_, int bottom_) : top{top_}, bottom{bottom_} {}
};

class TopBottomPeriodicBoundary {
 public:
  /**
   * @brief Construct a new TopBottomPeriodicBoundary object
   *
   * @param cells Cell index pairs of top and bottom boundaries
   */
  TopBottomPeriodicBoundary(const std::vector<TopBottomIndexPair>& cells)
      : cells_{cells} {}

  /**
   * @brief Construct a new Top Bottom Periodic Boundary object
   *
   * @tparam F
   * @param f Functor to calculate cell index pairs of top and bottom boundaries
   */
  template <typename F>
  TopBottomPeriodicBoundary(F&& f) : cells_{f()} {}

  /**
   * @brief Apply boundary condition
   *
   * @tparam T
   * @param f Distribution function
   */
  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    for (auto [top, bottom] : top_) {
      // clang-format off
      // f(1, top) = f(1, bottom);
      f(2, bottom) = f(2, top);
      // f(3, top) = f(3, bottom);
      f(4, top)    = f(4, bottom);
      f(5, bottom) = f(5, top);
      f(6, bottom) = f(6, top);
      f(7, top)    = f(7, bottom);
      f(8, top)    = f(8, bottom);
      // clang-format on
    }
  }

 private:
  /// A list of cell index pairs of top & bottom boundaries.
  std::vector<TopBottomIndexPair> cells_;
};

struct LeftRightIndexPair {
  int left;   ///< Cell index of left boundary
  int right;  ///< Cell index of right boundary

  LeftRightIndexPair(int left_, int right_) : left{left_}, right{right_} {}
};

class LeftRightPeriodicBoundary {
 public:
  /**
   * @brief Construct a new LeftRightPeriodicBoundary object
   *
   * @param cells Cell indices of left and right boundaries
   */
  LeftRightPeriodicBoundary(const std::vector<LeftRightIndexPair>& cells)
      : cells_{cells} {}

  /**
   * @brief Construct a new LeftRightPeriodicBoundary object
   *
   * @tparam F
   * @param f Functor to calculate cell index pairs of left and right
   * boundaries.
   */
  template <typename F>
  LeftRightPeriodicBoundary(F&& f) : cells_{f()} {}

  /**
   * @brief Apply boundary condition.
   *
   * @tparam T
   * @param f Distribution function
   */
  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    for (auto [left, right] : cells_) {
      // clang-format off
      f(1, left)  = f(1, right);
      // f(2, left) = fold(2, right);
      f(3, right) = f(3, left);
      // f(4, left) = fold(4, right);
      f(5, left)  = f(5, right);
      f(6, right) = f(6, left);
      f(7, right) = f(7, left);
      f(8, left)  = f(8, right);
      // clang-format on
    }
  }

 private:
  /// A list of cell index pairs of left & right boundaries.
  std::vector<LeftRightIndexPair> cells_;
};

}  // namespace lbm

#endif  // LBM_PERIODIC_BOUNDARY_HPP