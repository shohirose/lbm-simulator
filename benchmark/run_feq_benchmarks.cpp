#include <benchmark/benchmark.h>

#include <Eigen/Core>

class ComponentMajorLayoutFixture : public ::benchmark::Fixture {
 public:
  void SetUp(const ::benchmark::State& st) {
    // clang-format off
    c_ << 0,  1,  0, -1,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1;
    w_ << 4.0 / 9.0,
          1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
    // clang-format on
    const auto n = st.range(0);
    feq_ = Eigen::Matrix<double, 9, Eigen::Dynamic>::Zero(9, n);
    u_ = Eigen::Matrix<double, 2, Eigen::Dynamic>::Zero(2, n);
    rho_ = Eigen::VectorXd::Ones(n);
  }

  void TearDown(const ::benchmark::State& state) {}

  Eigen::Matrix<double, 2, 9> c_;
  Eigen::Matrix<double, 9, 1> w_;
  Eigen::Matrix<double, 2, Eigen::Dynamic> u_;
  Eigen::VectorXd rho_;
  Eigen::Matrix<double, 9, Eigen::Dynamic> feq_;
};

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod1Test)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu = c_.transpose() * u_;
    const Eigen::Matrix<double, 1, Eigen::Dynamic> u2 =
        u_.colwise().squaredNorm();
    feq_ = ((w_ * rho_.transpose()).array() *
            (1.0 + 3.0 * cu.array() + 4.5 * cu.array().square() -
             1.5 * u2.replicate<9, 1>().array()))
               .matrix();
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod1Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod2Test)
(benchmark::State& st) {
  for (auto _ : st) {
    for (int i = 0; i < feq_.cols(); ++i) {
      const auto u2 = u_.col(i).squaredNorm();
      for (int k = 0; k < 9; ++k) {
        const auto cu = c_.col(k).dot(u_.col(i));
        feq_(k, i) =
            w_(k) * rho_(i) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod2Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod3Test)
(benchmark::State& st) {
  for (auto _ : st) {
    for (int i = 0; i < feq_.cols(); ++i) {
      const Eigen::Vector2d ui = u_.col(i);
      const auto u2 = ui.squaredNorm();
      for (int k = 0; k < 9; ++k) {
        const auto cu = c_.col(k).dot(ui);
        feq_(k, i) =
            w_(k) * rho_(i) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod3Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod4Test)
(benchmark::State& st) {
  for (auto _ : st) {
    for (int i = 0; i < feq_.cols(); ++i) {
      const auto u2 = u_.col(i).squaredNorm();
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u_.col(i);
      for (int k = 0; k < 9; ++k) {
        feq_(k, i) = w_(k) * rho_(i) *
                     (1.0 + 3.0 * cu(k) + 4.5 * cu(k) * cu(k) - 1.5 * u2);
      }
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod4Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod5Test)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    for (int i = 0; i < feq_.cols(); ++i) {
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u_.col(i);
      for (int k = 0; k < 9; ++k) {
        feq_(k, i) = w_(k) * rho_(i) *
                     (1.0 + 3.0 * cu(k) + 4.5 * cu(k) * cu(k) - 1.5 * u2(i));
      }
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod5Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod6Test)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    for (int i = 0; i < feq_.cols(); ++i) {
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u_.col(i);
      feq_.col(i) =
          (rho_(i) * w_.array() *
           (1.0 + 3.0 * cu.array() + 4.5 * cu.array().square() - 1.5 * u2(i)))
              .matrix();
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod6Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod7Test)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu = c_.transpose() * u_;
    for (int i = 0; i < feq_.cols(); ++i) {
      feq_.col(i) = (rho_(i) * w_.array() *
                     (1.0 + 3.0 * cu.col(i).array() +
                      4.5 * cu.col(i).array().square() - 1.5 * u2(i)))
                        .matrix();
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod7Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(ComponentMajorLayoutFixture, FeqMethod8Test)
(benchmark::State& st) {
  for (auto _ : st) {
    for (int i = 0; i < feq_.cols(); ++i) {
      const Eigen::Vector2d u = u_.col(i);
      const auto u2 = u.squaredNorm();
      feq_(0, i) = rho_(i) * w_(0) * (1.0 - 1.5 * u2);
      feq_(1, i) = rho_(i) * w_(1) *
                   (1.0 + 3.0 * u.x() + 4.5 * u.x() * u.x() - 1.5 * u2);
      feq_(2, i) = rho_(i) * w_(2) *
                   (1.0 + 3.0 * u.y() + 4.5 * u.y() * u.y() - 1.5 * u2);
      feq_(3, i) = rho_(i) * w_(3) *
                   (1.0 - 3.0 * u.x() + 4.5 * u.x() * u.x() - 1.5 * u2);
      feq_(4, i) = rho_(i) * w_(4) *
                   (1.0 - 3.0 * u.y() + 4.5 * u.y() * u.y() - 1.5 * u2);
      {
        const auto cu = u.x() + u.y();
        feq_(5, i) =
            rho_(i) * w_(5) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
      {
        const auto cu = -u.x() + u.y();
        feq_(6, i) =
            rho_(i) * w_(6) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
      {
        const auto cu = -u.x() - u.y();
        feq_(7, i) =
            rho_(i) * w_(7) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
      {
        const auto cu = u.x() - u.y();
        feq_(8, i) =
            rho_(i) * w_(8) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
      }
    }
  }
}

BENCHMARK_REGISTER_F(ComponentMajorLayoutFixture, FeqMethod8Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

// class CellMajorLayoutFixture : public ::benchmark::Fixture {
//  public:
//   void SetUp(const ::benchmark::State& st) {
//     // clang-format off
//     c_.transpose() << 0,  1,  0, -1,  0,  1, -1, -1,  1,
//                       0,  0,  1,  0, -1,  1,  1, -1, -1;
//     w_ << 4.0 / 9.0,
//           1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
//           1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
//     // clang-format on
//     const auto n = st.range(0);
//     feq_ = Eigen::Matrix<double, Eigen::Dynamic, 9>::Zero(n, 9);
//     u_ = Eigen::Matrix<double, Eigen::Dynamic, 2>::Zero(n, 2);
//     rho_ = Eigen::VectorXd::Ones(n);
//   }

//   void TearDown(const ::benchmark::State& st) {}

//   Eigen::Matrix<double, 9, 2> c_;
//   Eigen::Matrix<double, 9, 1> w_;
//   Eigen::Matrix<double, Eigen::Dynamic, 2> u_;
//   Eigen::VectorXd rho_;
//   Eigen::Matrix<double, Eigen::Dynamic, 9> feq_;
// };

// BENCHMARK_DEFINE_F(CellMajorLayoutFixture, FeqMethod1Test)
// (benchmark::State& st) {
//   for (auto _ : st) {
//     for (int k = 0; k < 9; ++k) {
//       for (int i = 0; i < feq_.rows(); ++i) {
//         const auto u2 = u_.row(i).squaredNorm();
//         const auto cu = c_.row(k).dot(u_.row(i));
//         feq_(i, k) =
//             w_(k) * rho_(i) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
//       }
//     }
//   }
// }

// BENCHMARK_DEFINE_F(CellMajorLayoutFixture, FeqMethod2Test)
// (benchmark::State& st) {
//   for (auto _ : st) {
//     Eigen::Map<const Eigen::VectorXd> ux(&u_(0, 0), u_.rows());
//     Eigen::Map<const Eigen::VectorXd> uy(&u_(0, 1), u_.rows());
//     Eigen::Map<const Eigen::Matrix<double, 9, 1>> cx(&c_(0, 0));
//     Eigen::Map<const Eigen::Matrix<double, 9, 1>> cy(&c_(0, 1));
//     for (int k = 0; k < 9; ++k) {
//       for (int i = 0; i < feq_.rows(); ++i) {
//         const auto u2 = ux(i) * ux(i) + uy(i) * uy(i);
//         const auto cu = cx(k) * ux(i) + cy(k) * uy(i);
//         feq_(i, k) =
//             w_(k) * rho_(i) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
//       }
//     }
//   }
// }

// BENCHMARK_DEFINE_F(CellMajorLayoutFixture, FeqMethod3Test)
// (benchmark::State& st) {
//   for (auto _ : st) {
//     Eigen::Map<const Eigen::VectorXd> ux(&u_(0, 0), u_.rows());
//     Eigen::Map<const Eigen::VectorXd> uy(&u_(0, 1), u_.rows());
//     Eigen::Map<const Eigen::Matrix<double, 9, 1>> cx(&c_(0, 0));
//     Eigen::Map<const Eigen::Matrix<double, 9, 1>> cy(&c_(0, 1));
//     const Eigen::VectorXd u2 =
//         (ux.array().square() + uy.array().square()).matrix();
//     for (int k = 0; k < 9; ++k) {
//       const Eigen::VectorXd cu = cx(k) * ux + cy(k) * uy;
//       feq_.col(k) = (w_(k) * rho_.array() *
//                      (1.0 + 3.0 * cu.array() + 4.5 * cu.array().square() -
//                       1.5 * u2.array()))
//                         .matrix();
//     }
//   }
// }

// BENCHMARK_DEFINE_F(CellMajorLayoutFixture, FeqMethod4Test)
// (benchmark::State& st) {
//   for (auto _ : st) {
//     Eigen::Map<const Eigen::VectorXd> ux(&u_(0, 0), u_.rows());
//     Eigen::Map<const Eigen::VectorXd> uy(&u_(0, 1), u_.rows());
//     Eigen::Map<const Eigen::Matrix<double, 9, 1>> cx(&c_(0, 0));
//     Eigen::Map<const Eigen::Matrix<double, 9, 1>> cy(&c_(0, 1));
//     const Eigen::VectorXd u2 =
//         (ux.array().square() + uy.array().square()).matrix();
//     for (int k = 0; k < 9; ++k) {
//       const Eigen::VectorXd cu = cx(k) * ux + cy(k) * uy;
//       feq_.col(k) = (w_(k) * rho_.array() *
//                      (1.0 + 3.0 * (cx(k) * ux + cy(k) * uy).array() +
//                       4.5 * (cx(k) * ux + cy(k) * uy).array().square() -
//                       1.5 * u2.array()))
//                         .matrix();
//     }
//   }
// }

// BENCHMARK_REGISTER_F(CellMajorLayoutFixture, FeqMethod1Test)
//     ->RangeMultiplier(4)
//     ->Range(64, 256 << 8);

// BENCHMARK_REGISTER_F(CellMajorLayoutFixture, FeqMethod2Test)
//     ->RangeMultiplier(4)
//     ->Range(64, 256 << 8);

// BENCHMARK_REGISTER_F(CellMajorLayoutFixture, FeqMethod3Test)
//     ->RangeMultiplier(4)
//     ->Range(64, 256 << 8);

// BENCHMARK_REGISTER_F(CellMajorLayoutFixture, FeqMethod4Test)
//     ->RangeMultiplier(4)
//     ->Range(64, 256 << 8);

BENCHMARK_MAIN();