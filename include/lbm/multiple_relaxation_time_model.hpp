#ifndef LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP
#define LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP

#include <Eigen/Core>
#include <array>

namespace lbm {

struct MultipleRelaxationTimeModelParameters {
  std::array<double, 9> relaxation_parameters;

  MultipleRelaxationTimeModelParameters() = default;

  MultipleRelaxationTimeModelParameters(
      const std::array<double, 9>& relaxation_parameters_)
      : relaxation_parameters{relaxation_parameters_} {}

  double get_relaxation_time() const noexcept {
    return 1.0 / relaxation_parameters[7];
  }
};

class MultipleRelaxationTimeModel {
 public:
  using Matrix9d = Eigen::Matrix<double, 9, 9>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;

  MultipleRelaxationTimeModel(
      const MultipleRelaxationTimeModelParameters& params);

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
             [[maybe_unused]] const Eigen::MatrixBase<T3>& u) const noexcept {
    f -= C_ * (f - feq);
  }

 private:
  Matrix9d C_;  ///< = M^{-1} * S * M
};

}  // namespace lbm

#endif  // LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP