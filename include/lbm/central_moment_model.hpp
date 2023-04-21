#ifndef LBM_CENTRAL_MOMENT_MODEL_HPP
#define LBM_CENTRAL_MOMENT_MODEL_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <array>

namespace lbm {

struct CentralMomentModelParameters {
  std::array<double, 9> relaxation_matrix;

  CentralMomentModelParameters() = default;

  CentralMomentModelParameters(const std::array<double, 9>& s_)
      : relaxation_matrix{s_} {}

  CentralMomentModelParameters(const CentralMomentModelParameters&) = default;
  CentralMomentModelParameters(CentralMomentModelParameters&&) = default;

  CentralMomentModelParameters& operator=(
      const std::array<double, 9>& s) noexcept {
    relaxation_matrix = s;
  }
  CentralMomentModelParameters& operator=(const CentralMomentModelParameters&) =
      default;
  CentralMomentModelParameters& operator=(CentralMomentModelParameters&&) =
      default;

  auto* data() noexcept { return relaxation_matrix.data(); }
  const auto* data() const noexcept { return relaxation_matrix.data(); }

  auto& operator[](int i) noexcept { return relaxation_matrix[i]; }
  const auto& operator[](int i) const noexcept { return relaxation_matrix[i]; }
};

class CentralMomentModel {
 public:
  using Matrix9d = Eigen::Matrix<double, 9, 9>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;

  CentralMomentModel(const CentralMomentModelParameters& params);

  template <typename T1, typename T2, typename T3>
  void apply(Eigen::MatrixBase<T1>& f, const Eigen::MatrixBase<T2>& feq,
             const Eigen::MatrixBase<T3>& u) const noexcept {
    for (int i = 0; i < f.cols(); ++i) {
      const Matrix9d N = get_N(u(0, i), u(1, i));
      const Matrix9d Ninv = get_Ninv(u(0, i), u(1, i));
      Vector9d df = feq.col(i) - f.col(i);
      df = M_ * df;
      df = N * df;
      df = S_.asDiagonal() * df;
      df = Ninv * df;
      df = Minv_ * df;
      f.col(i) += df;
    }
  }

 private:
  static Matrix9d get_N(double ux, double uy) noexcept;

  static Matrix9d get_Ninv(double ux, double uy) noexcept;

  Matrix9d M_;
  Matrix9d Minv_;
  Vector9d S_;
};

}  // namespace lbm

#endif  // LBM_CENTRAL_MOMENT_MODEL_HPP