#ifndef LBM_CENTRAL_MOMENT_MODEL_HPP
#define LBM_CENTRAL_MOMENT_MODEL_HPP

#include <fmt/core.h>

#include <Eigen/Core>
#include <array>

namespace lbm {

struct CentralMomentModelParameters {
  std::array<double, 9> relaxation_parameters;

  CentralMomentModelParameters() = default;

  CentralMomentModelParameters(const std::array<double, 9>& s_)
      : relaxation_parameters{s_} {}

  double get_relaxation_time() const noexcept {
    return 1.0 / relaxation_parameters[4];
  }
};

class CentralMomentModel {
 public:
  using Matrix9d = Eigen::Matrix<double, 9, 9>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;

  CentralMomentModel(const CentralMomentModelParameters& params);

  /**
   * @brief Apply collision model
   *
   * @tparam T1
   * @tparam T2
   * @tparam T3
   * @param f Distribution function
   * @param feq Equilibrium distribution function
   * @param u Velocity
   */
  template <typename T1, typename T2, typename T3>
  void apply(Eigen::MatrixBase<T1>& f, const Eigen::MatrixBase<T2>& feq,
             const Eigen::MatrixBase<T3>& u) const noexcept {
    for (int i = 0; i < f.cols(); ++i) {
      const Matrix9d N = get_shift_matrix(u(0, i), u(1, i));
      const Matrix9d Ninv = get_inverse_shift_matrix(u(0, i), u(1, i));
      Vector9d df = feq.col(i) - f.col(i);
      // Transformation from distribution functions to raw moments
      Vector9d dm = M_ * df;
      // Transformation from raw moments to central moments
      Vector9d dmc = N * dm;
      // Relaxation
      dmc = S_.asDiagonal() * dmc;
      // Back transformation from central moments to raw moments
      dm = Ninv * dmc;
      // Back transformation from raw moments to distribution funcitons
      df = Minv_ * dm;
      // Update distribution functions
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