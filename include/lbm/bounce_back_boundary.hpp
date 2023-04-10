#ifndef LBM_BOUNCE_BACK_BOUNDARY_HPP
#define LBM_BOUNCE_BACK_BOUNDARY_HPP

#include <Eigen/Core>
#include <variant>
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

/**
 * @brief Cell index tuple
 *
 */
struct IndexTuple {
  int i;   ///< x direction
  int j;   ///< y direction
  int id;  ///< cell index

  /**
   * @brief Construct a new IndexTuple object
   *
   * @param i_ Index in the x direction
   * @param j_ Index in the y direction
   * @param id_ Cell index
   */
  IndexTuple(int i_, int j_, int id_) : i{i_}, j{j_}, id{id_} {}
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
        f(2, cell.id) = f(4, cell.id);
        f(5, cell.id) = f(7, cell.id);
        f(6, cell.id) = f(8, cell.id);
      }
    } else if constexpr (B == BoundaryNormal::Down) {
      for (auto&& cell : cells_) {
        f(4, cell.id) = f(2, cell.id);
        f(7, cell.id) = f(5, cell.id);
        f(8, cell.id) = f(6, cell.id);
      }
    } else if constexpr (B == BoundaryNormal::Left) {
      for (auto&& cell : cells_) {
        f(3, cell.id) = f(1, cell.id);
        f(6, cell.id) = f(8, cell.id);
        f(7, cell.id) = f(5, cell.id);
      }
    } else if constexpr (B == BoundaryNormal::Right) {
      for (auto&& cell : cells_) {
        f(1, cell.id) = f(3, cell.id);
        f(8, cell.id) = f(6, cell.id);
        f(5, cell.id) = f(7, cell.id);
      }
    }
  }

 private:
  std::vector<IndexTuple> cells_;  ///< A list of cells on the boundary.
};

}  // namespace lbm

#endif  // LBM_BOUNCE_BACK_BOUNDARY_HPP