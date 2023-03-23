#ifndef LBM_CARTESIAN_GRID_2D_HPP
#define LBM_CARTESIAN_GRID_2D_HPP

#include <array>
#include <cassert>

namespace lbm {

class CartesianGrid2d {
 public:
  CartesianGrid2d(int nx, int ny) : shape_{nx, ny} {}

  CartesianGrid2d(const std::array<int, 2>& shape) : shape_{shape} {}

  int size() const noexcept { return shape_[0] * shape_[1]; }

  int shape(int i) const noexcept {
    assert(i >= 0 && i < shape_.size() && "out of bounds error!");
    return shape_[i];
  }

  int nx() const noexcept {
    return this->shape(0);
  }

  int ny() const noexcept {
    return this->shape(1);
  }

  const auto& shape() const noexcept { return shape_; }

  int index(int i, int j) const noexcept {
    assert(i >= 0 && i <= shape_[1] && "out of bounds error!");
    assert(j >= 0 && j <= shape_[0] && "out of bounds error!");
    return i * shape_[0] + j;
  }

  int periodic_index(int i, int j) const noexcept {
    if (i < 0) {
      i += shape_[1];
    }
    if (j < 0) {
      j += shape_[0];
    }
    if (i >= shape_[1]) {
      i -= shape_[1];
    }
    if (j >= shape_[0]) {
      j -= shape_[0];
    }
    return this->index(i, j);
  }

 private:
  std::array<int, 2> shape_;
};

}  // namespace lbm

#endif  // LBM_CARTESIAN_GRID_2D_HPP