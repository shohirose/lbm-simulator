#include <benchmark/benchmark.h>

#include <Eigen/Core>

class EquilibriumDistributionFunctionFixture : public ::benchmark::Fixture {
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

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, NaiveImplTest)
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

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture, NaiveImplTest)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, ForLoopTest1)
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

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture, ForLoopTest1)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, ForLoopTest2)
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

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture, ForLoopTest2)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest1)
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

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest1)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest2)
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

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest2)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest3)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    for (int i = 0; i < feq_.cols(); ++i) {
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u_.col(i);
      const Eigen::Matrix<double, 9, 1> cu2 = cu.array().square().matrix();
      feq_.col(i) = (rho_(i) * w_.array() *
                     (1.0 + 3.0 * cu.array() + 4.5 * cu2.array() - 1.5 * u2(i)))
                        .matrix();
    }
  }
}

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest3)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest4)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    for (int i = 0; i < feq_.cols(); ++i) {
      const Eigen::Matrix<double, 9, 1> cu = c_.transpose() * u_.col(i);
      const Eigen::Matrix<double, 9, 1> cu2 = cu.array().square().matrix();
      for (int k = 0; k < 9; ++k) {
        feq_(k, i) =
            rho_(i) * w_(k) * (1.0 + 3.0 * cu(k) + 4.5 * cu2(k) - 1.5 * u2(i));
      }
    }
  }
}

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest4)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest5)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu = c_.transpose() * u_;
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu2 =
        cu.array().square().matrix();
    for (int i = 0; i < feq_.cols(); ++i) {
      for (int k = 0; k < 9; ++k) {
        feq_(k, i) = rho_(i) * w_(k) *
                     (1.0 + 3.0 * cu(k, i) + 4.5 * cu2(k, i) - 1.5 * u2(i));
      }
    }
  }
}

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest5)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest6)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu = c_.transpose() * u_;
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu2 =
        cu.array().square().matrix();
    for (int i = 0; i < feq_.cols(); ++i) {
      feq_.col(i) = (rho_(i) * w_.array() *
                     (1.0 + 3.0 * cu.col(i).array() + 4.5 * cu2.col(i).array() -
                      1.5 * u2(i)))
                        .matrix();
    }
  }
}

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest6)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(EquilibriumDistributionFunctionFixture, EigenVectorizeTest7)
(benchmark::State& st) {
  for (auto _ : st) {
    const Eigen::VectorXd u2 = u_.colwise().squaredNorm().transpose();
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu = c_.transpose() * u_;
    const Eigen::Matrix<double, 9, Eigen::Dynamic> cu2 =
        cu.array().square().matrix();
    feq_ =
        (rho_.replicate<9, 1>().array() * w_.replicate(1, feq_.cols()).array() *
         (1.0 + 3.0 * cu.array() + 4.5 * cu2.array() -
          1.5 * u2.replicate<9, 1>().array()))
            .matrix();
  }
}

BENCHMARK_REGISTER_F(EquilibriumDistributionFunctionFixture,
                     EigenVectorizeTest7)
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