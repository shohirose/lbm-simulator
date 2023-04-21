#ifndef LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP
#define LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

namespace lbm {

struct MultipleRelaxationTimeModelParameters {
  double se;
  double sq;
  double seps;
  double tau;

  MultipleRelaxationTimeModelParameters() = default;

  MultipleRelaxationTimeModelParameters(double se_, double sq_, double seps_,
                                        double tau_)
      : se{se_}, sq{sq_}, seps{seps_}, tau{tau_} {}
};

class MultipleRelaxationTimeModel {
 public:
  using Matrix9d = Eigen::Matrix<double, 9, 9>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;

  MultipleRelaxationTimeModel(
      const MultipleRelaxationTimeModelParameters& params);

  template <typename T1, typename T2>
  void apply(Eigen::MatrixBase<T1>& f,
             const Eigen::MatrixBase<T2>& feq) const noexcept {
    f -= C_ * (f - feq);
  }

 private:
  Matrix9d C_;  ///< = M^{-1} * S * M
};

}  // namespace lbm

#endif  // LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP