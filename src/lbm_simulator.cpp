#include <CLI/CLI.hpp>

#include "lbm/cavity_flow_simulator.hpp"
#include "lbm/poiseuille_flow_simulator.hpp"

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
      lbm::PoiseuilleFlowSimulator simulator(filename1);
      simulator.run();
    } else if (*sub2) {
      lbm::CavityFlowSimulator simulator(filename2);
      simulator.run();
    }
  } catch (const std::exception& e) {
    fmt::print(stderr, "Error: {}", e.what());
    std::exit(EXIT_FAILURE);
  }
}