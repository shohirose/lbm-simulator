#ifndef LBM_PROPAGATOR_HPP
#define LBM_PROPAGATOR_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

/**
 * @brief Propagator of distribution function for all cells.
 *
 */
class AllCellPropagator {
 public:
  AllCellPropagator(const CartesianGrid2d& grid) : grid_{grid} {}

  template <typename T>
  void apply(Eigen::MatrixBase<T>& f,
             Eigen::MatrixBase<T>& fold) const noexcept {
    fold = f;
    for (int j = 0; j < grid_.ny(); ++j) {
      for (int i = 0; i < grid_.nx(); ++i) {
        const auto n = grid_.index(i, j);
        if (i + 1 < grid_.nx()) {
          f(1, grid_.index(i + 1, j)) = fold(1, n);
        }
        if (j + 1 < grid_.ny()) {
          f(2, grid_.index(i, j + 1)) = fold(2, n);
        }
        if (i - 1 >= 0) {
          f(3, grid_.index(i - 1, j)) = fold(3, n);
        }
        if (j - 1 >= 0) {
          f(4, grid_.index(i, j - 1)) = fold(4, n);
        }
        if (i + 1 < grid_.nx() && j + 1 < grid_.ny()) {
          f(5, grid_.index(i + 1, j + 1)) = fold(5, n);
        }
        if (i - 1 >= 0 && j + 1 < grid_.ny()) {
          f(6, grid_.index(i - 1, j + 1)) = fold(6, n);
        }
        if (i - 1 >= 0 && j - 1 >= 0) {
          f(7, grid_.index(i - 1, j - 1)) = fold(7, n);
        }
        if (i + 1 < grid_.nx() && j - 1 >= 0) {
          f(8, grid_.index(i + 1, j - 1)) = fold(8, n);
        }
      }
    }
  }

 private:
  CartesianGrid2d grid_;
};

/**
 * @brief Propagator of distribution function for internal cells.
 *
 */
class InternalCellPropagator {
 public:
  InternalCellPropagator(const CartesianGrid2d& grid) : grid_{grid} {}

  template <typename T>
  void apply(Eigen::MatrixBase<T>& f,
             Eigen::MatrixBase<T>& fold) const noexcept {
    fold = f;
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f(1, grid_.index(i + 1, j    )) = fold(1, n);
        f(2, grid_.index(i,     j + 1)) = fold(2, n);
        f(3, grid_.index(i - 1, j    )) = fold(3, n);
        f(4, grid_.index(i    , j - 1)) = fold(4, n);
        f(5, grid_.index(i + 1, j + 1)) = fold(5, n);
        f(6, grid_.index(i - 1, j + 1)) = fold(6, n);
        f(7, grid_.index(i - 1, j - 1)) = fold(7, n);
        f(8, grid_.index(i + 1, j - 1)) = fold(8, n);
        // clang-format on
      }
    }
  }

 private:
  CartesianGrid2d grid_;
};

}  // namespace lbm

#endif  // LBM_PROPAGATOR_HPP