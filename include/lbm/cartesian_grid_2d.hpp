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
  CartesianGrid2d(int ni, int nj) : shape_{ni, nj} {}

  CartesianGrid2d(const std::array<int, 2>& shape) : shape_{shape} {}

  int size() const noexcept { return shape_[0] * shape_[1]; }

  int shape(int i) const noexcept {
    assert(i >= 0 && i < shape_.size() && "out of bounds error!");
    return shape_[i];
  }

  int ni() const noexcept { return this->shape(0); }

  int nj() const noexcept { return this->shape(1); }

  const auto& shape() const noexcept { return shape_; }

  int index(int i, int j) const noexcept {
    assert(i >= 0 && i <= shape_[0] && "out of bounds error!");
    assert(j >= 0 && j <= shape_[1] && "out of bounds error!");
    return i * shape_[1] + j;
  }

  template <typename Boundary,
            std::enable_if_t<std::is_same_v<Boundary, TopBoundary> ||
                                 std::is_same_v<Boundary, BottomBoundary>,
                             std::nullptr_t> = nullptr>
  int index(Boundary boundary, int j) const noexcept {
    return detail::boundary_index_impl<Boundary>::eval(j, shape_[0], shape_[1]);
  }

  template <typename Boundary,
            std::enable_if_t<std::is_same_v<Boundary, LeftBoundary> ||
                                 std::is_same_v<Boundary, RightBoundary>,
                             std::nullptr_t> = nullptr>
  int index(int i, Boundary boundary) const noexcept {
    return detail::boundary_index_impl<Boundary>::eval(i, shape_[0], shape_[1]);
  }

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