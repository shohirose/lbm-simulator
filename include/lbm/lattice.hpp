#ifndef LBM_LATTICE_HPP
#define LBM_LATTICE_HPP

#include <Eigen/Core>

namespace lbm {

enum class LatticeType {
  D2Q5,
  D2Q9,
  D3Q15,
  D3Q19,
  D3Q27,
};

template <LatticeType T>
struct Lattice {};

template <>
struct Lattice<LatticeType::D2Q5> {
  static constexpr auto dimensions = 2;
  static constexpr auto connections = 5;

  using LatticeVector = Eigen::Matrix<double, dimensions, connections>;
  using WeightVector = Eigen::Matrix<double, connections, 1>;

  static LatticeVector get_lattice_vector() noexcept {
    LatticeVector c;
    // clang-format off
    c << 0,  1,  0, -1,  0,
         0,  0,  1,  0, -1;
    // clang-format on
    return c;
  }

  static WeightVector get_weight() noexcept {
    WeightVector w;
    // clang-format off
    w << 4.0 / 9.0, 1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0;
    // clang-format on
    return w;
  }
};

template <>
struct Lattice<LatticeType::D2Q9> {
  static constexpr auto dimensions = 2;
  static constexpr auto connections = 9;

  using LatticeVector = Eigen::Matrix<double, dimensions, connections>;
  using WeightVector = Eigen::Matrix<double, connections, 1>;

  static LatticeVector get_lattice_vector() noexcept {
    LatticeVector c;
    // clang-format off
    c << 0,  1,  0, -1,  0,  1, -1, -1,  1,
         0,  0,  1,  0, -1,  1,  1, -1, -1;
    // clang-format on
    return c;
  }

  static WeightVector get_weight() noexcept {
    WeightVector w;
    // clang-format off
    w << 4.0 / 9.0,
         1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
         1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0;
    // clang-format on
    return w;
  }
};

template <>
struct Lattice<LatticeType::D3Q15> {
  static constexpr auto dimensions = 3;
  static constexpr auto connections = 15;

  using LatticeVector = Eigen::Matrix<double, dimensions, connections>;
  using WeightVector = Eigen::Matrix<double, connections, 1>;

  static LatticeVector get_lattice_vector() noexcept {
    LatticeVector c;
    // clang-format off
    c << 0,  1,  0,  0, -1,  0,  0,  1, -1,  1,  1,  1, -1, -1, -1,
         0,  0,  1,  0,  0, -1,  0,  1,  1, -1, -1, -1, -1,  1, -1,
         0,  0,  0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1;
    // clang-format on
    return c;
  }

  static WeightVector get_weight() noexcept {
    WeightVector w;
    // clang-format off
    w << 2.0 / 9.0,
         1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
         1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0,
         1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0;
    // clang-format on
    return w;
  }
};

template <>
struct Lattice<LatticeType::D3Q19> {
  static constexpr auto dimensions = 3;
  static constexpr auto connections = 19;

  using LatticeVector = Eigen::Matrix<double, dimensions, connections>;
  using WeightVector = Eigen::Matrix<double, connections, 1>;

  static LatticeVector get_lattice_vector() noexcept {
    LatticeVector c;
    // clang-format off
    c.transpose()
      << 0,  0,  0,
         1,  0,  0,
        -1,  0,  0,
         0,  1,  0,
         0, -1,  0,
         0,  0,  1,
         0,  0, -1,
         1,  1,  0,
        -1, -1,  0,
         1, -1,  0,
        -1,  1,  0,
         1,  0,  1,
        -1,  0, -1,
         1,  0, -1,
        -1,  0,  1,
         0,  1,  1,
         0, -1, -1,
         0,  1, -1,
         0, -1,  1;
    // clang-format on
    return c;
  }

  static WeightVector get_weight() noexcept {
    WeightVector w;
    // clang-format off
    w << 1.0 / 3.0,
         1.0 / 18.0,
         1.0 / 18.0,
         1.0 / 18.0,
         1.0 / 18.0,
         1.0 / 18.0,
         1.0 / 18.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0,
         1.0 / 36.0;
    // clang-format on
    return w;
  }
};

template <>
struct Lattice<LatticeType::D3Q27> {
  static constexpr auto dimensions = 3;
  static constexpr auto connections = 27;

  using LatticeVector = Eigen::Matrix<double, dimensions, connections>;
  using WeightVector = Eigen::Matrix<double, connections, 1>;

  static LatticeVector get_lattice_vector() noexcept {
    LatticeVector c;
    // clang-format off
    c.transpose()
      <<  0,  0,  0,
          1,  0,  0,
         -1,  0,  0,
          0,  1,  0,
          0, -1,  0,
          0,  0,  1,
          0,  0, -1,
          1,  1,  0,
         -1, -1,  0,
          1, -1,  0,
         -1,  1,  0,
          1,  0,  1,
         -1,  0, -1,
          1,  0, -1,
         -1,  0,  1,
          0,  1,  1,
          0, -1, -1,
          0,  1, -1,
          0, -1,  1,
          1,  1,  1,
         -1, -1, -1,
         -1,  1,  1,
          1, -1, -1,
          1, -1,  1,
         -1,  1, -1,
          1,  1, -1,
         -1, -1,  1;
    // clang-format on
    return c;
  }

  static WeightVector get_weight() noexcept {
    WeightVector w;
    // clang-format off
    w << 8.0 / 27.0,
         2.0 / 27.0,
         2.0 / 27.0,
         2.0 / 27.0,
         2.0 / 27.0,
         2.0 / 27.0,
         2.0 / 27.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 54.0,
         1.0 / 216.0,
         1.0 / 216.0,
         1.0 / 216.0,
         1.0 / 216.0,
         1.0 / 216.0,
         1.0 / 216.0,
         1.0 / 216.0,
         1.0 / 216.0;
    // clang-format on
    return w;
  }
};

}  // namespace lbm

#endif  // LBM_LATTICE_HPP