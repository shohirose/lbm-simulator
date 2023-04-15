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
      : C_{} {
    Eigen::Matrix<double, 9, 9> M;
    Eigen::Matrix<double, 9, 9> Minv;
    Eigen::Matrix<double, 9, 1> S;
    // clang-format off
    M <<  1,  1,  1,  1,  1,  1,  1,  1,  1,
         -4, -1, -1, -1, -1,  2,  2,  2,  2,
          4, -2, -2, -2, -2,  1,  1,  1,  1,
          0,  1,  0, -1,  0,  1, -1, -1,  1,
          0, -2,  0,  2,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1,
          0,  0, -2,  0,  2,  1,  1, -1, -1,
          0,  1, -1,  1, -1,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  1, -1,  1, -1;
    // clang-format on
    Minv = M.inverse();
    const auto snu = 1 / params.tau;
    S << 0, params.se, params.seps, 0, params.sq, 0, params.sq, snu, snu;
    C_ = Minv * S.asDiagonal() * M;
  }

  template <typename T1, typename T2>
  void apply(Eigen::MatrixBase<T1>& f,
             const Eigen::MatrixBase<T2>& feq) const noexcept {
    f -= C_ * (f - feq);
  }

 private:
  Eigen::Matrix<double, 9, 9> C_;  ///< = M^{-1} * S * M
};

}  // namespace lbm

#endif  // LBM_MULTIPLE_RELAXATION_TIME_MODEL_HPP