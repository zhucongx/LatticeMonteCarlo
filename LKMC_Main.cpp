#include <chrono>
#include "Home.h"
#include "StateChangePredictor.h"
int main(int argc, char *argv[]) {
  api::Parameter parameter;
  if (argc == 1) {
    std::cout << "No input parameter filename. Opening lkmc_param.txt" << std::endl;
    parameter = api::Parameter("lkmc_param.txt");
  } else {
    parameter = api::Parameter(argc, argv);
  }

  parameter.PrintParameters();

  if (parameter.method == "FirstKmc") {
    auto first_kmc = api::BuildFirstKmcMpiFromParameter(parameter);
    first_kmc.Simulate();
  } else if (parameter.method == "ChainKmc") {
    auto chain_kmc = api::BuildChainKmcMpiFromParameter(parameter);
    chain_kmc.Simulate();
  } else if (parameter.method == "Cluster") {

  } else if (parameter.method == "SimulatedAnnealing") {
    auto simulated_annealing_mc = api::BuildSimulatedAnnealingFromParameter(parameter);
    simulated_annealing_mc.Simulate();
  }

  // ansys::Iterator test(0, 1e4,
  //                      Element("Al"),
  //                      {Element("Al"),
  //                       Element("Mg"),
  //                       Element("Zn")},
  //                      4, "quartic_coefficients.json");
  // // test.SerialRunReformat();
  // test.SerialRunCluster();



//   pred::TotalEnergyPredictor a("quartic_coefficients.json",
//                                std::set<Element>{Element("Al"), Element("Mg"),
//                                                  Element("Zn")});
//   const auto conf0 = cfg::GenerateFCC({10, 10, 10}, Element("Al"));
//   size_t Zn, Mg;
//   std::cout << "Mg Zn Energy" << std::endl;
// #pragma omp parallel for default(none) shared(conf0, a, std::cout) private(Zn, Mg)
//   for (Mg = 0; Mg <= 200; ++Mg) {
//     for (Zn = 0; Zn <= 200; ++Zn) {
//       if (Mg + Zn > 200) continue;
//       if (Mg + Zn < 4) continue;
//       auto conf1 = GenerateSoluteConfigFromExcitingPure(
//           conf0,
//           {{Element("Mg"), Mg},
//            {Element("Zn"), Zn}});
//       double energy = a.GetEnergy(conf1);
//       std::string output = std::to_string(Mg) + ' ' +
//           std::to_string(Zn) + ' ' + std::to_string(energy) + '\n';
//       std::cout << output << std::flush;
//     }
//   }
//
//   return 0;
}
