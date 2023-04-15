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
  MultipleRelaxationTimeModel(
      const MultipleRelaxationTimeModelParameters& params)
      : M_{}, Minv_{}, S_{}, C_{} {
    // clang-format off
    M_ <<  1,  1,  1,  1,  1,  1,  1,  1,  1,
          -4, -1, -1, -1, -1,  2,  2,  2,  2,
           4, -2, -2, -2, -2,  1,  1,  1,  1,
           0,  1,  0, -1,  0,  1, -1, -1,  1,
           0, -2,  0,  2,  0,  1, -1, -1,  1,
           0,  0,  1,  0, -1,  1,  1, -1, -1,
           0,  0, -2,  0,  2,  1,  1, -1, -1,
           0,  1, -1,  1, -1,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  1, -1,  1, -1;
    // clang-format on
    Minv_ = M_.inverse();
    const auto snu = 1 / params.tau;
    S_ << 0, params.se, params.seps, 0, params.sq, 0, params.sq, snu, snu;
    C_ = Minv_ * S_.asDiagonal() * M_;
  }

  template <typename T1, typename T2>
  void apply(Eigen::MatrixBase<T1>& f,
             const Eigen::MatrixBase<T2>& feq) const noexcept {
    f -= C_ * (f - feq);
  }

 private:
  Eigen::Matrix<double, 9, 9> M_;     ///< Transformation matrix
  Eigen::Matrix<double, 9, 9> Minv_;  ///< Inverse of M
  Eigen::Matrix<double, 9, 1> S_;     ///< Relaxation time vector
  Eigen::Matrix<double, 9, 9> C_;     ///< = Minv_ * S_.asDiagonal() * M_
};

}  // namespace lbm

#endif  // LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP