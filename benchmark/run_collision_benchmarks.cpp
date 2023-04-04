#include <benchmark/benchmark.h>

#include <Eigen/Core>

class CollisionProcessFixture : public ::benchmark::Fixture {
 public:
  void SetUp(const ::benchmark::State& st) {
    const auto n = st.range(0);
    const auto size = n * n;
    f_ = Eigen::Matrix<double, 9, Eigen::Dynamic>::Ones(9, size);
    feq_ = Eigen::Matrix<double, 9, Eigen::Dynamic>::Ones(9, size);
    tau_ = 0.52;
  }
  void TearDown(const ::benchmark::State& st) {}

  Eigen::Matrix<double, 9, 1> f_;
  Eigen::Matrix<double, 9, 1> feq_;
  double tau_;
};

BENCHMARK_DEFINE_F(CollisionProcessFixture, CollisionMethod1Test)
(benchmark::State& st) {
  for (auto _ : st) {
    f_ = f_ - (f_ - feq_) / tau_;
  }
}

BENCHMARK_REGISTER_F(CollisionProcessFixture, CollisionMethod1Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(CollisionProcessFixture, CollisionMethod2Test)
(benchmark::State& st) {
  for (auto _ : st) {
    f_ *= 1.0 - 1.0 / tau_;
    f_ += feq_ / tau_;
  }
}

BENCHMARK_REGISTER_F(CollisionProcessFixture, CollisionMethod2Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(CollisionProcessFixture, CollisionMethod3Test)
(benchmark::State& st) {
  for (auto _ : st) {
    f_ *= 1.0 - 1.0 / tau_;
    f_.noalias() += feq_ / tau_;
  }
}

BENCHMARK_REGISTER_F(CollisionProcessFixture, CollisionMethod3Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);

BENCHMARK_DEFINE_F(CollisionProcessFixture, CollisionMethod4Test)
(benchmark::State& st) {
  for (auto _ : st) {
    f_ *= 1.0 - 1.0 / tau_;
    f_.noalias() += feq_ * (1.0 / tau_);
  }
}

BENCHMARK_REGISTER_F(CollisionProcessFixture, CollisionMethod4Test)
    ->RangeMultiplier(4)
    ->Range(64, 256 << 8);