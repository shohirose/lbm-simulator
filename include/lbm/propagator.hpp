#ifndef LBM_PROPAGATOR_HPP
#define LBM_PROPAGATOR_HPP

#include <Eigen/Core>
#include <vector>

#include "lbm/cartesian_grid_2d.hpp"

namespace lbm {

struct PropagationIndexTuple {
  int lattice_id;
  int cell_from;
  int cell_to;

  PropagationIndexTuple(int id_, int from_, int to_)
      : lattice_id{id_}, cell_from{from_}, cell_to{to_} {}
};

class Propagator {
 public:
  Propagator(const CartesianGrid2d& grid) : connections_{} {
    connections_.reserve(grid.size() * 9);
    for (int j = 0; j < grid.ny(); ++j) {
      for (int i = 0; i < grid.nx(); ++i) {
        const auto from = grid.index(i, j);
        if (i + 1 < grid.nx()) {
          connections_.emplace_back(1, from, grid.index(i + 1, j));
        }
        if (j + 1 < grid.ny()) {
          connections_.emplace_back(2, from, grid.index(i, j + 1));
        }
        if (i - 1 >= 0) {
          connections_.emplace_back(3, from, grid.index(i - 1, j));
        }
        if (j - 1 >= 0) {
          connections_.emplace_back(4, from, grid.index(i, j - 1));
        }
        if (i + 1 < grid.nx() && j + 1 < grid.ny()) {
          connections_.emplace_back(5, from, grid.index(i + 1, j + 1));
        }
        if (i - 1 >= 0 && j + 1 < grid.ny()) {
          connections_.emplace_back(6, from, grid.index(i - 1, j + 1));
        }
        if (i - 1 >= 0 && j - 1 >= 0) {
          connections_.emplace_back(7, from, grid.index(i - 1, j - 1));
        }
        if (i + 1 < grid.nx() && j - 1 >= 0) {
          connections_.emplace_back(8, from, grid.index(i + 1, j - 1));
        }
      }
    }
  }

  template <typename T>
  void apply(Eigen::MatrixBase<T>& f,
             Eigen::MatrixBase<T>& fold) const noexcept {
    fold = f;
    for (auto&& [k, from, to] : connections_) {
      f(k, to) = fold(k, from);
    }
  }

 private:
  std::vector<PropagationIndexTuple> connections_;
};

}  // namespace lbm

#endif  // LBM_PROPAGATOR_HPP