#ifndef LBM_CARTESIAN_GRID_2D_HPP
#define LBM_CARTESIAN_GRID_2D_HPP

#include <array>
#include <cassert>

namespace lbm {

class TopBoundary {};
class BottomBoundary {};
class LeftBoundary {};
class RightBoundary {};

static constexpr auto top = TopBoundary{};
static constexpr auto bottom = BottomBoundary{};
static constexpr auto left = LeftBoundary{};
static constexpr auto right = RightBoundary{};

namespace detail {

template <typename Boundary>
struct boundary_index_impl {};

template <>
struct boundary_index_impl<TopBoundary> {
  static int eval(int j, int ni, int nj) noexcept { return (ni - 1) * nj + j; }
};

template <>
struct boundary_index_impl<BottomBoundary> {
  static int eval(int j, int ni, int nj) noexcept { return j; }
};

template <>
struct boundary_index_impl<LeftBoundary> {
  static int eval(int i, int ni, int nj) noexcept { return i * nj; }
};

template <>
struct boundary_index_impl<RightBoundary> {
  static int eval(int i, int ni, int nj) noexcept { return i * nj + nj - 1; }
};

}  // namespace detail

class CartesianGrid2d {
 public:
  /**
   * @brief Construct a new Cartesian Grid 2d object
   *
   * @param ni Number of cells in I direction
   * @param nj Number of cells in J direction
   */
  CartesianGrid2d(int ni, int nj) : shape_{ni, nj} {}

  /**
   * @brief Construct a new Cartesian Grid 2d object
   *
   * @param shape Shape of grid: (ni, nj)
   */
  CartesianGrid2d(const std::array<int, 2>& shape) : shape_{shape} {}

  /// @brief Return the total number of cells
  int size() const noexcept { return shape_[0] * shape_[1]; }

  /**
   * @brief Return the number of cells in the given direction
   *
   * @param i Direction. I = 0, J = 1.
   */
  int shape(int i) const noexcept {
    assert(i >= 0 && i < shape_.size() && "out of bounds error!");
    return shape_[i];
  }

  int ni() const noexcept { return this->shape(0); }

  int nj() const noexcept { return this->shape(1); }

  const auto& shape() const noexcept { return shape_; }

  /**
   * @brief Return the cell index of a given I and J.
   *
   * @param i Index in the I direction
   * @param j Index in the J direction
   */
  int index(int i, int j) const noexcept {
    assert(i >= 0 && i <= shape_[0] && "out of bounds error!");
    assert(j >= 0 && j <= shape_[1] && "out of bounds error!");
    return i * shape_[1] + j;
  }

  /**
   * @brief Return the cell index along either top or bottom boundaries.
   *
   * @tparam Boundary TopBoundary or BottomBoundary
   * @param boundary lbm::top or lbm::bottom
   * @param j Index in the J direction
   */
  template <typename Boundary,
            std::enable_if_t<std::is_same_v<Boundary, TopBoundary> ||
                                 std::is_same_v<Boundary, BottomBoundary>,
                             std::nullptr_t> = nullptr>
  int index(Boundary boundary, int j) const noexcept {
    return detail::boundary_index_impl<Boundary>::eval(j, shape_[0], shape_[1]);
  }

  /**
   * @brief Return the cell index along either left or right boundaries.
   *
   * @tparam Boundary LeftBoundary or RightBoundary
   * @param i Index in the I direction
   * @param boundary lbm::left or lbm::right
   */
  template <typename Boundary,
            std::enable_if_t<std::is_same_v<Boundary, LeftBoundary> ||
                                 std::is_same_v<Boundary, RightBoundary>,
                             std::nullptr_t> = nullptr>
  int index(int i, Boundary boundary) const noexcept {
    return detail::boundary_index_impl<Boundary>::eval(i, shape_[0], shape_[1]);
  }

  /**
   * @brief Return the periodic cell index.
   *
   * If i or j exceeds its bounds, either ni or nj is added or subtracted to fit
   * in the bounds.
   *
   * @param i Index in the I direction
   * @param j Index in the J direction
   */
  int periodic_index(int i, int j) const noexcept {
    if (i < 0) {
      i += shape_[0];
    }
    if (j < 0) {
      j += shape_[1];
    }
    if (i >= shape_[0]) {
      i -= shape_[0];
    }
    if (j >= shape_[1]) {
      j -= shape_[1];
    }
    return this->index(i, j);
  }

 private:
  std::array<int, 2> shape_;
};

}  // namespace lbm

#endif  // LBM_CARTESIAN_GRID_2D_HPP