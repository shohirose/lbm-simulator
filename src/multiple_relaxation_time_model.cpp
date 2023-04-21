#include "lbm\multiple_relaxation_time_model.hpp"

namespace lbm {

MultipleRelaxationTimeModel::MultipleRelaxationTimeModel(
    const MultipleRelaxationTimeModelParameters& params)
    : C_{} {
  Matrix9d M;
  Matrix9d Minv;
  Vector9d S;
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

}  // namespace lbm
