#ifndef LBM_CARTESIAN_GRID_2D_HPP
#define LBM_CARTESIAN_GRID_2D_HPP

#include <array>
#include <cassert>

namespace lbm {

class CartesianGrid2d {
 public:
  CartesianGrid2d() = default;
  
  /**
   * @brief Construct a new Cartesian Grid 2d object
   *
   * @param nx Number of cells in the x direction
   * @param ny Number of cells in the y direction
   */
  CartesianGrid2d(int nx, int ny) : shape_{nx, ny} {}

  /**
   * @brief Construct a new Cartesian Grid 2d object
   *
   * @param shape Shape of grid: (nx, ny)
   */
  CartesianGrid2d(const std::array<int, 2>& shape) : shape_{shape} {}

  /// @brief Return the total number of cells
  int size() const noexcept { return shape_[0] * shape_[1]; }

  /**
   * @brief Return the number of cells in the given direction
   *
   * @param i Direction. X = 0, Y = 1.
   */
  int shape(int i) const noexcept {
    assert(i >= 0 && i < shape_.size() && "out of bounds error!");
    return shape_[i];
  }

  int nx() const noexcept { return this->shape(0); }

  int ny() const noexcept { return this->shape(1); }

  const auto& shape() const noexcept { return shape_; }

  /**
   * @brief Return the cell index.
   *
   * @param i Index in the x direction
   * @param j Index in the y direction
   */
  int index(int i, int j) const noexcept {
    assert(i >= 0 && i < shape_[0] && "Error: out of bounds!");
    assert(j >= 0 && j < shape_[1] && "Error: out of bounds!");
    return i + j * shape_[0];
  }

  /**
   * @brief Return the periodic cell index.
   *
   * If i or j exceeds its bounds, either ni or nj is added or subtracted to fit
   * in the bounds.
   *
   * @param i Index in the x direction
   * @param j Index in the y direction
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