#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include "lbm/cavity_flow_simulator.hpp"
#include "lbm/poiseuille_flow_simulator.hpp"

using Json = nlohmann::json;
namespace fs = std::filesystem;

namespace lbm {

/**
 * @brief Conversion function from json to CavityFlowParameter
 */
void from_json(const Json& j, CavityFlowParameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("wallVelocity").get_to(params.wall_velocity);
  j.at("reynoldsNumber").get_to(params.reynolds_number);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
  j.at("outputDirectory").get_to(params.output_directory);
}

/**
 * @brief Conversion function from json to PoiseuilleFlowParameter
 */
void from_json(const Json& j, PoiseuilleFlowParameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("externalForce").get_to(params.external_force);
  j.at("relaxationTime").get_to(params.relaxation_time);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
  j.at("outputDirectory").get_to(params.output_directory);
}

}  // namespace lbm

struct PoiseuilleProblem {
  using Parameters = lbm::PoiseuilleFlowParameters;
  using Simulator = lbm::PoiseuilleFlowSimulator;
};

struct CavityProblem {
  using Parameters = lbm::CavityFlowParameters;
  using Simulator = lbm::CavityFlowSimulator;
};

template <typename Problem>
void run_simulator(const std::string& filename) {
  using Parameters = typename Problem::Parameters;
  using Simulator = typename Problem::Simulator;

  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error(
        fmt::format("Error: could not open a file: {}", filename));
  }
  const auto j = Json::parse(file);
  const auto params = j.get<Parameters>();
  Simulator simulator(params);
  simulator.run();
}

int main(int argc, char* argv[]) {
  CLI::App app{"Lattice Boltzaman Simulator"};

  auto* sub1 =
      app.add_subcommand("poiseuille", "Run 2-D Poiseuille flow simulator.");
  std::string filename1 = "poiseuille.json";
  sub1->add_option("-f,--file", filename1,
                   "Input JSON file. (default = 'poiseuille.json')");

  auto* sub2 = app.add_subcommand("cavity", "Run 2-D Cavity flow simulator.");
  std::string filename2 = "cavity.json";
  sub2->add_option("-f,--file", filename2,
                   "Input JSON file. (default = 'cavity.json')");

  CLI11_PARSE(app, argc, argv);

  try {
    if (*sub1) {
      run_simulator<PoiseuilleProblem>(filename1);
    } else if (*sub2) {
      run_simulator<CavityProblem>(filename2);
    }
  } catch (const std::exception& e) {
    fmt::print(stderr, "Error: {}", e.what());
    std::exit(EXIT_FAILURE);
  }
}