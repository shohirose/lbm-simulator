#include <CLI/CLI.hpp>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "lbm/poiseuille_flow_simulator.hpp"

namespace fs = std::filesystem;

using Eigen::VectorXd;
using json = nlohmann::json;
using Parameters = lbm::PoiseuilleFlowSimulator::Parameters;

namespace lbm {

void from_json(const json& j, Parameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("externalForce").get_to(params.external_force);
  j.at("relaxationTime").get_to(params.relaxation_time);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
}

}  // namespace lbm

void check_path(const fs::path& p);
Parameters get_parameters(const fs::path& p);

int main(int argc, char* argv[]) {
  CLI::App app{"2D Poiseuille flow simulator"};
  std::string filename = "poiseuille-flow-2d.json";
  app.add_option("-f,--file", filename, "Input file name in JSON format.");

  std::string dir = "result/poiseuille-flow-2d";
  app.add_option("-o,--output-dir", dir, "Output directory path.");

  CLI11_PARSE(app, argc, argv);

  fs::path p(filename);
  check_path(p);
  const auto params = get_parameters(p);

  lbm::PoiseuilleFlowSimulator simulator(params);
  const VectorXd ux = simulator.calc_velocity();

  fs::path pd(dir);
  try {
    create_directories(pd);
    std::ofstream file(pd / fs::path("ux.txt"));
    file << ux << std::endl;
  } catch (const fs::filesystem_error& e) {
    fmt::print(stderr, "Error occured while creating directories: {}\n",
               e.what());
    std::exit(EXIT_FAILURE);
  }
}

void check_path(const fs::path& p) {
  if (!fs::exists(p)) {
    fmt::print(stderr, "Error: file does not exists: {}\n", p.string());
    std::exit(EXIT_FAILURE);
  }
  if (!fs::is_regular_file(p)) {
    fmt::print(stderr, "Error: not a regulalr file: {}\n", p.string());
    std::exit(EXIT_FAILURE);
  }
  if (p.extension() != ".json") {
    fmt::print(stderr, "Error: not a JSON file: {}\n", p.string());
    std::exit(EXIT_FAILURE);
  }
}

Parameters get_parameters(const fs::path& p) {
  try {
    std::ifstream file(p.string());
    const auto j = json::parse(file);
    return j.get<Parameters>();
  } catch (std::exception& e) {
    fmt::print(stderr, "Error occured while reading a JSON file: {}\n  {}",
               p.string(), e.what());
    std::exit(EXIT_FAILURE);
  }
}