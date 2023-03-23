#include <filesystem>
#include <fstream>
#include <iostream>

#include "lbm/poiseuille_flow_simulator.hpp"

namespace fs = std::filesystem;

using Eigen::VectorXd;

int main() {
  const std::array<int, 2> shape = {21, 20};
  const std::array<double, 2> external_force = {0.00001, 0.0};
  const double relaxation_time = 0.56;
  const double error_limit = 1e-10;
  const int print_freq = 5000;
  const int max_iter = 1000000;
  const lbm::Parameters params{shape,       external_force, relaxation_time,
                               error_limit, print_freq,     max_iter};
  lbm::PoiseuilleFlowSimulator simulator(params);
  const VectorXd ux = simulator.calc_velocity();

  std::ofstream file(fs::path("ux.txt"));
  file << ux << std::endl;
}