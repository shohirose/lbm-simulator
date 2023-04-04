#include <benchmark/benchmark.h>

#include <Eigen/Core>

#include "lbm/cartesian_grid_2d.hpp"

class PropagationProcessFixture : public ::benchmark::Fixture {
 public:
  void SetUp(const ::benchmark::State& st) {
    const auto n = st.range(0);
    grid_ = lbm::CartesianGrid2d(n, n);
    f_ = Eigen::Matrix<double, 9, Eigen::Dynamic>::Zero(9, grid_.size());
    fold_ = Eigen::Matrix<double, 9, Eigen::Dynamic>::Zero(9, grid_.size());
  }

  void TearDown(const ::benchmark::State& st) {}

  lbm::CartesianGrid2d grid_;
  Eigen::Matrix<double, 9, Eigen::Dynamic> f_;
  Eigen::Matrix<double, 9, Eigen::Dynamic> fold_;
};

BENCHMARK_DEFINE_F(PropagationProcessFixture, PropagationMethod1Test)
(benchmark::State& st) {
  for (auto _ : st) {
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f_(1, grid_.index(i + 1, j    )) = fold_(1, n);
        f_(2, grid_.index(i    , j + 1)) = fold_(2, n);
        f_(3, grid_.index(i - 1, j    )) = fold_(3, n);
        f_(4, grid_.index(i    , j - 1)) = fold_(4, n);
        f_(5, grid_.index(i + 1, j + 1)) = fold_(5, n);
        f_(6, grid_.index(i - 1, j + 1)) = fold_(6, n);
        f_(7, grid_.index(i - 1, j - 1)) = fold_(7, n);
        f_(8, grid_.index(i + 1, j - 1)) = fold_(8, n);
        // clang-format on
      }
    }
  }
}

BENCHMARK_DEFINE_F(PropagationProcessFixture, PropagationMethod2Test)
(benchmark::State& st) {
  for (auto _ : st) {
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f_(1, grid_.index(i + 1, j    )) = fold_(1, n);
        f_(3, grid_.index(i - 1, j    )) = fold_(3, n);
        // clang-format on
      }
    }
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f_(1, grid_.index(i + 1, j    )) = fold_(1, n);
        f_(5, grid_.index(i + 1, j + 1)) = fold_(5, n);
        f_(8, grid_.index(i + 1, j - 1)) = fold_(8, n);
        // clang-format on
      }
    }
    for (int j = 1; j < grid_.ny() - 1; ++j) {
      for (int i = 1; i < grid_.nx() - 1; ++i) {
        const auto n = grid_.index(i, j);
        // clang-format off
        f_(3, grid_.index(i - 1, j    )) = fold_(3, n);
        f_(6, grid_.index(i - 1, j + 1)) = fold_(6, n);
        f_(7, grid_.index(i - 1, j - 1)) = fold_(7, n);
        // clang-format on
      }
    }
  }
}

BENCHMARK_REGISTER_F(PropagationProcessFixture, PropagationMethod1Test)
    ->RangeMultiplier(2)
    ->Range(16, 16 << 4);

BENCHMARK_REGISTER_F(PropagationProcessFixture, PropagationMethod2Test)
    ->RangeMultiplier(2)
    ->Range(16, 16 << 4);