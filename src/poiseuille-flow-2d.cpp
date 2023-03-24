#include <CLI/CLI.hpp>
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

void to_json(json& j, const Parameters& params) {
  j = json{{"gridShape", params.grid_shape},
           {"externalForce", params.external_force},
           {"relaxationTime", params.relaxation_time},
           {"errorLimit", params.error_limit},
           {"printFrequency", params.print_frequency},
           {"maxIteration", params.max_iter}};
}

void from_json(const json& j, Parameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("externalForce").get_to(params.external_force);
  j.at("relaxationTime").get_to(params.relaxation_time);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
}

}  // namespace lbm

int main(int argc, char* argv[]) {
  CLI::App app{"2D Poiseuille flow simulator"};
  std::string filename = "poiseuille-flow-2d.json";
  app.add_option("-f,--file", filename, "Input file name in JSON format.");

  CLI11_PARSE(app, argc, argv);

  fs::path p(filename);
  if (!fs::exists(p)) {
    std::cerr << fmt::format("Error: file does not exists: {}\n", p.string());
    std::exit(EXIT_FAILURE);
  }
  if (!fs::is_regular_file(p)) {
    std::cerr << fmt::format("Error: not a regulalr file: {}\n", p.string());
    std::exit(EXIT_FAILURE);
  }
  if (p.extension() != ".json") {
    std::cerr << fmt::format("Error: not a JSON file: {}\n", p.string());
    std::exit(EXIT_FAILURE);
  }

  json j;
  {
    std::ifstream file(p.string());
    try {
      j = json::parse(file);
    } catch (std::exception& e) {
      std::cerr << fmt::format(
          "Error occured while reading a JSON file: {}\n  {}", p.string(),
          e.what());
      std::exit(EXIT_FAILURE);
    }
  }
  const auto params = j.get<Parameters>();

  lbm::PoiseuilleFlowSimulator simulator(params);
  const VectorXd ux = simulator.calc_velocity();

  std::ofstream file(fs::path("ux.txt"));
  file << ux << std::endl;
}