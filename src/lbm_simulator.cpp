#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include "lbm/cavity_flow_simulator.hpp"
#include "lbm/poiseuille_flow_simulator.hpp"

using Json = nlohmann::json;
namespace fs = std::filesystem;

namespace lbm {

void from_json(const Json& j, SingleRelaxationTimeModelParameters& params) {
  j.at("tau").get_to(params.tau);
}

void from_json(const Json& j, MultipleRelaxationTimeModelParameters& params) {
  j.at("se").get_to(params.se);
  j.at("sq").get_to(params.sq);
  j.at("seps").get_to(params.seps);
  j.at("tau").get_to(params.tau);
}

void from_json(const Json& j, CentralMomentModelParameters& params) {
  j.at("s").get_to(params.s);
}

/**
 * @brief Conversion function from json to CavityFlowParameter
 */
void from_json(const Json& j, CavityFlowParameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("wallVelocity").get_to(params.wall_velocity);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
  j.at("outputDirectory").get_to(params.output_directory);
  if (j.contains("singleRelaxationTimeModel")) {
    params.collision_params = j.at("singleRelaxationTimeModel")
                                  .get<SingleRelaxationTimeModelParameters>();
  } else if (j.contains("multipleRelaxationTimeModel")) {
    params.collision_params = j.at("multipleRelaxationTimeModel")
                                  .get<MultipleRelaxationTimeModelParameters>();
  } else if (j.contains("centralMomentModel")) {
    params.collision_params =
        j.at("centralMomentModel").get<CentralMomentModelParameters>();
  } else {
    throw std::runtime_error(
        "Error: collision model parameters not found: "
        "[singleRelaxationTimeModel, multipleRelaxationTimeModel, "
        "centralMomentModel]");
  }
}

/**
 * @brief Conversion function from json to PoiseuilleFlowParameter
 */
void from_json(const Json& j, PoiseuilleFlowParameters& params) {
  j.at("gridShape").get_to(params.grid_shape);
  j.at("externalForce").get_to(params.external_force);
  j.at("errorLimit").get_to(params.error_limit);
  j.at("printFrequency").get_to(params.print_frequency);
  j.at("maxIteration").get_to(params.max_iter);
  j.at("outputDirectory").get_to(params.output_directory);
  if (j.contains("singleRelaxationTimeModel")) {
    params.collision_params = j.at("singleRelaxationTimeModel")
                                  .get<SingleRelaxationTimeModelParameters>();
  } else if (j.contains("multipleRelaxationTimeModel")) {
    params.collision_params = j.at("multipleRelaxationTimeModel")
                                  .get<MultipleRelaxationTimeModelParameters>();
  } else if (j.contains("centralMomentModel")) {
    params.collision_params =
        j.at("centralMomentModel").get<CentralMomentModelParameters>();
  } else {
    throw std::runtime_error(
        "Error: collision model parameters not found: "
        "[singleRelaxationTimeModel, multipleRelaxationTimeModel, "
        "centralMomentModel]");
  }
}

}  // namespace lbm

auto get_input_data(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error(
        fmt::format("Error: could not open a file: {}", filename));
  }
  return Json::parse(file);
}

template <typename Simulator>
void run_simulator(const Json& j) {
  const auto params = j.get<typename Simulator::Parameters>();
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
      const auto j = get_input_data(filename1);
      run_simulator<lbm::PoiseuilleFlowSimulator>(j);
    } else if (*sub2) {
      const auto j = get_input_data(filename2);
      run_simulator<lbm::CavityFlowSimulator>(j);
    } else {
      throw std::runtime_error(
          "Error: subcommand is missing: [poiseuille, cavity]");
    }
  } catch (const std::exception& e) {
    fmt::print(stderr, "Error: {}", e.what());
    std::exit(EXIT_FAILURE);
  }
}