#include <filesystem>
#include <fstream>
#include <iostream>

#include "lbm/poiseuille_flow_simulator.hpp"

namespace fs = std::filesystem;

using Eigen::MatrixXd, Eigen::Map, Eigen::Stride, Eigen::Dynamic,
    Eigen::Unaligned;

int main() {
  const std::array<int, 2> shape = {21, 21};
  const std::array<double, 2> external_force = {0.00001, 0.0};
  const double relaxation_time = 0.56;
  const double error_limit = 1e-10;
  const int print_freq = 5000;
  const int max_iter = 1000000;
  const lbm::Params params{shape,       external_force, relaxation_time,
                           error_limit, print_freq,     max_iter};
  lbm::PoiseuilleFlowSimulator simulator(params);
  const MatrixXd u = simulator.calc_velocity();

  {
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> ux(
        &u(0, 0), shape[0], shape[1], Stride<Dynamic, 2>(shape[1]*2, 2));

    std::ofstream file(fs::path("ux.txt"));
    file << ux.transpose() << std::endl;
  }

  {
    Map<const MatrixXd, Unaligned, Stride<Dynamic, 2>> uy(
        &u(1, 0), shape[0], shape[1], Stride<Dynamic, 2>(shape[1]*2, 2));

    std::ofstream file(fs::path("uy.txt"));
    file << uy.transpose() << std::endl;
  }
}