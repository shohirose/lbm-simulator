#ifndef LBM_NON_ORTHOGONAL_CENTRAL_MOMENT_MODEL_HPP
#define LBM_NON_ORTHOGONAL_CENTRAL_MOMENT_MODEL_HPP

#include <Eigen/Core>
#include <array>
#include <iostream>

namespace lbm {

struct NonOrthogonalCentralMomentModelParameters {
  std::array<double, 9> s;  ///< Relaxation matrix parameters

  NonOrthogonalCentralMomentModelParameters() = default;

  NonOrthogonalCentralMomentModelParameters(const std::array<double, 9>& s_)
      : s{s_} {}
};

class NonOrthogonalCentralMomentModel {
 public:
  using Matrix9d = Eigen::Matrix<double, 9, 9>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;

  NonOrthogonalCentralMomentModel(
      const NonOrthogonalCentralMomentModelParameters& params)
      : C_{} {
    Eigen::Map<const Vector9d> S(params.s.data());
    Matrix9d M;
    Matrix9d Minv;
    // clang-format off
    M  <<  1,  1,  1,  1,  1,  1,  1,  1,  1,
           0,  1,  0, -1,  0,  1, -1, -1,  1,
           0,  0,  1,  0, -1,  1,  1, -1, -1,
           0,  1,  1,  1,  1,  2,  2,  2,  2,
           0,  1, -1,  1, -1,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  1, -1,  1, -1,
           0,  0,  0,  0,  0,  1,  1, -1, -1,
           0,  0,  0,  0,  0,  1, -1, -1,  1,
           0,  0,  0,  0,  0,  1,  1,  1,  1;
    Minv  << 1,    0,    0,   -1,     0,     0,     0,     0,     1,
             0,  0.5,    0, 0.25,  0.25,     0,     0,  -0.5,  -0.5,
             0,    0,  0.5, 0.25, -0.25,     0,  -0.5,     0,  -0.5,
             0, -0.5,    0, 0.25,  0.25,     0,     0,   0.5,  -0.5,
             0,    0, -0.5, 0.25, -0.25,     0,   0.5,     0,  -0.5,
             0,    0,    0,    0,     0,  0.25,  0.25,  0.25,  0.25,
             0,    0,    0,    0,     0, -0.25,  0.25, -0.25,  0.25,
             0,    0,    0,    0,     0,  0.25, -0.25, -0.25,  0.25,
             0,    0,    0,    0,     0, -0.25, -0.25,  0.25,  0.25;
    // clang-format on
    C_ = Minv * S.asDiagonal() * M;
  }

  template <typename T1, typename T2, typename T3>
  void apply(Eigen::MatrixBase<T1>& f, const Eigen::MatrixBase<T2>& feq,
             const Eigen::MatrixBase<T3>& u) const {
    for (int i = 0; i < f.cols(); ++i) {
      const Matrix9d N = get_N(u(0, i), u(1, i));
      const Matrix9d Ninv = get_Ninv(u(0, i), u(1, i));
      const Matrix9d C = Ninv * C_ * N;
      f.col(i) -= C * (f.col(i) - feq.col(i));
    }
  }

 private:
  static Eigen::Matrix<double, 9, 9> get_N(double ux, double uy) {
    const auto ux2 = ux * ux;
    const auto uy2 = uy * uy;
    Matrix9d N = Matrix9d::Identity();
    N(1, 0) = -ux;
    N(2, 0) = uy;
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

  static Eigen::Matrix<double, 9, 9> get_Ninv(double ux, double uy) {
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

  Eigen::Matrix<double, 9, 9> C_;
};

}  // namespace lbm

#endif  // LBM_NON_ORTHOGONAL_CENTRAL_MOMENT_MODEL_HPP