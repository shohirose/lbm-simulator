#include "lbm\multiple_relaxation_time_model.hpp"

namespace lbm {

MultipleRelaxationTimeModel::MultipleRelaxationTimeModel(
    const MultipleRelaxationTimeModelParameters& params)
    : C_{} {
  Matrix9d M;
  Matrix9d Minv;
  Eigen::Map<const Vector9d> S(params.relaxation_parameters.data());
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
  Minv << 1.00, -1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
          1.00, -0.25, -0.50,  1.50, -1.50,  0.00,  0.00,  2.25,  0.00,
          1.00, -0.25, -0.50,  0.00,  0.00,  1.50, -1.50, -2.25,  0.00,
          1.00, -0.25, -0.50, -1.50,  1.50,  0.00,  0.00,  2.25,  0.00,
          1.00, -0.25, -0.50,  0.00,  0.00, -1.50,  1.50, -2.25,  0.00,
          1.00,  0.50,  0.25,  1.50,  0.75,  1.50,  0.75,  0.00,  2.25,
          1.00,  0.50,  0.25, -1.50, -0.75,  1.50,  0.75,  0.00, -2.25,
          1.00,  0.50,  0.25, -1.50, -0.75, -1.50, -0.75,  0.00,  2.25,
          1.00,  0.50,  0.25,  1.50,  0.75, -1.50, -0.75,  0.00, -2.25;
  Minv /= 9.0;
  // clang-format on
  C_ = Minv * S.asDiagonal() * M;
}

}  // namespace lbm