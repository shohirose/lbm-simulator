#ifndef LBM_SINGLE_RELAXATION_TIME_MODEL_HPP
#define LBM_SINGLE_RELAXATION_TIME_MODEL_HPP

#include <Eigen/Core>

namespace lbm {

struct SingleRelaxationTimeModelParameters {
  double tau;  ///< Relaxation time

  SingleRelaxationTimeModelParameters() = default;
  
  SingleRelaxationTimeModelParameters(double tau_) : tau{tau_} {}
};

class SingleRelaxationTimeModel {
 public:
  /**
   * @brief Construct a new SingleRelaxationTimeModel object
   *
   * @param relaxation_time Relaxation time
   */
  SingleRelaxationTimeModel(const SingleRelaxationTimeModelParameters& params)
      : tau_{params.tau} {}

  template <typename T1, typename T2>
  void apply(Eigen::MatrixBase<T1>& f,
             const Eigen::MatrixBase<T2>& feq) const noexcept {
    // f = f - (f - feq)/tau_;
    const auto alpha = 1 / tau_;
    f *= 1 - alpha;
    f.noalias() += feq * alpha;
  }

 private:
  double tau_;  ///< Relaxation time
};

}  // namespace lbm

#endif  // LBM_SINGLE_RELAXATION_TIME_MODEL_HPP