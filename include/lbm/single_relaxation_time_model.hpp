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
    f = f - (1 / tau_) * (f - feq);
  }

 private:
  double tau_;  ///< Relaxation time
};

}  // namespace lbm

#endif  // LBM_SINGLE_RELAXATION_TIME_MODEL_HPP