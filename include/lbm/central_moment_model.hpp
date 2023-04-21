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
      const Matrix9d N = get_shift_matrix(u(0, i), u(1, i));
      const Matrix9d Ninv = get_inverse_shift_matrix(u(0, i), u(1, i));
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
  /**
   * @brief Compute shift matrix from raw moments to central moments
   *
   * @param ux Velocity in the x direction
   * @param uy Velocity in the y direction
   * @return Matrix9d Shift matrix
   */
  static Matrix9d get_shift_matrix(double ux, double uy) noexcept;

  /**
   * @brief Compute the inverse of shift matrix
   *
   * @param ux Velocity in the x direction
   * @param uy Velocity in the y direction
   * @return Matrix9d Inverse of shift matrix
   */
  static Matrix9d get_inverse_shift_matrix(double ux, double uy) noexcept;

  /// @brief Transformation matrix from distribution functions to raw moments
  Matrix9d M_;
  /// @brief Inverse matrix of M
  Matrix9d Minv_;
  /// @brief Relaxation matrix
  Vector9d S_;
};

}  // namespace lbm

#endif  // LBM_CENTRAL_MOMENT_MODEL_HPP