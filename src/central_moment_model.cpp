#include "lbm/central_moment_model.hpp"

#include <fmt/core.h>

namespace lbm {

CentralMomentModel::CentralMomentModel(
    const CentralMomentModelParameters& params)
    : M_{}, Minv_{}, S_{} {
  Eigen::Map<const Vector9d> S(params.relaxation_parameters.data());
  S_ = S;
  if (S_(0) != 0 || S_(1) != 0 || S_(2) != 0) {
    S_(0) = 0.0;
    S_(1) = 0.0;
    S_(2) = 0.0;
    fmt::print(
        "Warning: s0, s1, and s2 in the relaxation matrix are set to zero!");
  }
  // clang-format off
  M_ <<   1,  1,  1,  1,  1,  1,  1,  1,  1,
          0,  1,  0, -1,  0,  1, -1, -1,  1,
          0,  0,  1,  0, -1,  1,  1, -1, -1,
          0,  1,  1,  1,  1,  2,  2,  2,  2,
          0,  1, -1,  1, -1,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  1, -1,  1, -1,
          0,  0,  0,  0,  0,  1,  1, -1, -1,
          0,  0,  0,  0,  0,  1, -1, -1,  1,
          0,  0,  0,  0,  0,  1,  1,  1,  1;
  Minv_ <<  1,    0,    0,   -1,     0,     0,     0,     0,     1,
            0,  0.5,    0, 0.25,  0.25,     0,     0,  -0.5,  -0.5,
            0,    0,  0.5, 0.25, -0.25,     0,  -0.5,     0,  -0.5,
            0, -0.5,    0, 0.25,  0.25,     0,     0,   0.5,  -0.5,
            0,    0, -0.5, 0.25, -0.25,     0,   0.5,     0,  -0.5,
            0,    0,    0,    0,     0,  0.25,  0.25,  0.25,  0.25,
            0,    0,    0,    0,     0, -0.25,  0.25, -0.25,  0.25,
            0,    0,    0,    0,     0,  0.25, -0.25, -0.25,  0.25,
            0,    0,    0,    0,     0, -0.25, -0.25,  0.25,  0.25;
  // clang-format on
}

CentralMomentModel::Matrix9d CentralMomentModel::get_shift_matrix(
    double ux, double uy) noexcept {
  const auto ux2 = ux * ux;
  const auto uy2 = uy * uy;
  Matrix9d N = Matrix9d::Identity();
  N(1, 0) = -ux;
  N(2, 0) = -uy;
  N(3, 0) = ux2 + uy2;
  N(4, 0) = ux2 - uy2;
  N(5, 0) = ux * uy;
  N(6, 0) = -ux2 * uy;
  N(7, 0) = -ux * uy2;
  N(8, 0) = ux2 * uy2;
  N(3, 1) = -2 * ux;
  N(4, 1) = -2 * ux;
  N(5, 1) = -uy;
  N(6, 1) = 2 * ux * uy;
  N(7, 1) = uy2;
  N(8, 1) = -2 * ux * uy2;
  N(3, 2) = -2 * uy;
  N(4, 2) = 2 * uy;
  N(5, 2) = -ux;
  N(6, 2) = ux2;
  N(7, 2) = 2 * ux * uy;
  N(8, 2) = -2 * ux2 * uy;
  N(6, 3) = -0.5 * uy;
  N(7, 3) = -0.5 * ux;
  N(8, 3) = 0.5 * (ux2 + uy2);
  N(6, 4) = -0.5 * uy;
  N(7, 4) = 0.5 * ux;
  N(8, 4) = 0.5 * (-ux2 + uy2);
  N(6, 5) = -2 * ux;
  N(7, 5) = -2 * uy;
  N(8, 5) = 4 * ux * uy;
  N(8, 6) = -2 * uy;
  N(8, 7) = -2 * ux;
  return N;
}

CentralMomentModel::Matrix9d CentralMomentModel::get_inverse_shift_matrix(
    double ux, double uy) noexcept {
  const auto ux2 = ux * ux;
  const auto uy2 = uy * uy;
  Matrix9d Ninv = Matrix9d::Identity();
  Ninv(1, 0) = ux;
  Ninv(2, 0) = uy;
  Ninv(3, 0) = ux2 + uy2;
  Ninv(4, 0) = ux2 - uy2;
  Ninv(5, 0) = ux * uy;
  Ninv(6, 0) = ux2 * uy;
  Ninv(7, 0) = ux * uy2;
  Ninv(8, 0) = ux2 * uy2;
  Ninv(3, 1) = 2 * ux;
  Ninv(4, 1) = 2 * ux;
  Ninv(5, 1) = uy;
  Ninv(6, 1) = 2 * ux * uy;
  Ninv(7, 1) = uy2;
  Ninv(8, 1) = 2 * ux * uy2;
  Ninv(3, 2) = 2 * uy;
  Ninv(4, 2) = -2 * uy;
  Ninv(5, 2) = ux;
  Ninv(6, 2) = ux2;
  Ninv(7, 2) = 2 * ux * uy;
  Ninv(8, 2) = 2 * ux2 * uy;
  Ninv(6, 3) = 0.5 * uy;
  Ninv(7, 3) = 0.5 * ux;
  Ninv(8, 3) = 0.5 * (ux2 + uy2);
  Ninv(6, 4) = 0.5 * uy;
  Ninv(7, 4) = -0.5 * ux;
  Ninv(8, 4) = 0.5 * (-ux2 + uy2);
  Ninv(6, 5) = 2 * ux;
  Ninv(7, 5) = 2 * uy;
  Ninv(8, 5) = 4 * ux * uy;
  Ninv(8, 6) = 2 * uy;
  Ninv(8, 7) = 2 * ux;
  return Ninv;
}

}  // namespace lbm
