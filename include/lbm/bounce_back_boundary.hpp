#ifndef LBM_BOUNCE_BACK_BOUNDARY_HPP
#define LBM_BOUNCE_BACK_BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

namespace lbm {

/**
 * @brief Direction of boundary normal.
 *
 */
enum class BoundaryNormal {
  Up,
  Down,
  Left,
  Right,
};

template <BoundaryNormal B>
class BounceBackBoundary {
 public:
  BounceBackBoundary() = default;

  /**
   * @brief Construct a new BounceBackBoundary object
   *
   * @tparam F
   * @param f Functor creating a list of cells on the boundary
   */
  template <typename F>
  BounceBackBoundary(F&& f) : cells_{f()} {}

  /**
   * @brief Apply the boundary condition
   *
   * @tparam T
   * @param f Distribution function
   */
  template <typename T>
  void apply(Eigen::MatrixBase<T>& f) const noexcept {
    if constexpr (B == BoundaryNormal::Up) {
      for (auto&& cell : cells_) {
        f(2, cell) = f(4, cell);
        f(5, cell) = f(7, cell);
        f(6, cell) = f(8, cell);
      }
    } else if constexpr (B == BoundaryNormal::Down) {
      for (auto&& cell : cells_) {
        f(4, cell) = f(2, cell);
        f(7, cell) = f(5, cell);
        f(8, cell) = f(6, cell);
      }
    } else if constexpr (B == BoundaryNormal::Left) {
      for (auto&& cell : cells_) {
        f(3, cell) = f(1, cell);
        f(6, cell) = f(8, cell);
        f(7, cell) = f(5, cell);
      }
    } else if constexpr (B == BoundaryNormal::Right) {
      for (auto&& cell : cells_) {
        f(1, cell) = f(3, cell);
        f(8, cell) = f(6, cell);
        f(5, cell) = f(7, cell);
      }
    }
  }

 private:
  std::vector<int> cells_;  ///< A list of cells on the boundary.
};

}  // namespace lbm

#endif  // LBM_BOUNCE_BACK_BOUNDARY_HPP