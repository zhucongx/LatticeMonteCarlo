#include <chrono>
#include "Home.h"
#include "EnergyChangePredictorFaster.h"
int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  parameter.PrintParameters();

  if (parameter.method == "FirstKmc") {
    auto first_kmc = api::BuildFirstKmcMpiFromParameter(parameter);
    first_kmc.Simulate();
  } else if (parameter.method == "ChainKmcOmp") {
    auto chain_kmc_omp = api::BuildChainKmcOmpFromParameter(parameter);
    chain_kmc_omp.Simulate();
  } else if (parameter.method == "ChainKmcMpi") {
    auto chain_kmc_mpi = api::BuildChainKmcMpiFromParameter(parameter);
    chain_kmc_mpi.Simulate();
  } else if (parameter.method == "FindCluster") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunCluster();
  } else if (parameter.method == "Reformat") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunReformat();
  } else if (parameter.method == "SimulatedAnnealing") {
    auto simulated_annealing = api::BuildSimulatedAnnealingFromParameter(parameter);
    simulated_annealing.Simulate();
  } else if (parameter.method == "CanonicalMc") {
    auto canonical_mc = api::BuildCanonicalMcFromParameter(parameter);
    canonical_mc.Simulate();
  } else if (parameter.method == "CanonicalMcStepT") {
    auto canonical_mc_step_t = api::BuildCanonicalMcStepTFromParameter(parameter);
    canonical_mc_step_t.Simulate();
  } else {
    std::cout << "No such method: " << parameter.method << std::endl;
    return 1;
  }

  // pred::TotalEnergyPredictor a("quartic_coefficients.json",
  //                              std::set<Element>{Element("Al"), Element("Mg"),
  //                                                Element("Zn")});
  // auto conf = cfg::Config::ReadCfg("Large30_Mg25_Zn56.cfg");
  // double energy = a.GetEnergy(conf);
  // std::cout << energy << std::flush;

  // const auto conf0 = cfg::GenerateFCC({10, 10, 10}, Element("Al"));
  // size_t Zn, Mg;
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
  return 0;
}
