#include <chrono>
#include "Home.h"
#include "StateChangePredictor.h"
int main(int argc, char *argv[]) {
  // api::Parameter parameter;
  // if (argc == 1) {
  //   std::cout << "No input parameter filename. Opening kmc_param.txt" << std::endl;
  //   parameter = api::Parameter("kmc_param.txt");
  // } else {
  //   parameter = api::Parameter(argc, argv);
  // }
  // parameter.PrintParameters();
  // if (parameter.method == "First") {
  //   auto kmc = api::BuildFirstKmcMpiFromParameter(parameter);
  //   kmc.Simulate();
  //
  // } else if (parameter.method == "Chain") {
  //   auto kmc = api::BuildChainKmcMpiFromParameter(parameter);
  //   kmc.Simulate();
  // } else if (parameter.method == "Cluster") {
  //
  // } else if (parameter.method == "Metropolis") {
  //
  // }

//   ansys::Iterator test(0, 1e4,
//                        Element("Al"),
//                        {Element("Al"),
//                         Element("Mg"),
//                         Element("Zn")},
//                        4, "quartic_coefficients.json");
//   // test.SerialRunReformat();
//   test.SerialRunCluster();

  // auto t1 = std::chrono::high_resolution_clock::now();
  // pred::StateChangePredictor a("quartic_coefficients.json",
  //                              cfg::GenerateFCC(
  //                                  4.046, {8, 8, 8}, Element("Al")),
  //                              {Element("Al"),
  //                               Element("Mg"),
  //                               Element("Zn")});
  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << std::setprecision(16)
  //           << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << '\n';
  auto
  a = ansys::SimulatedAnnealing({8, 8, 8},
                                Element("Al"),
                                {{Element("Mg"), 20}, {Element("Zn"), 20},
                                 1e2, 1e5, 1e10, 1e6, "quartic_coefficients.json");
                                a.Simulate();


//   pred::TotalEnergyPredictor a("quartic_coefficients.json",
//                                std::set<Element>{Element("Al"), Element("Mg"),
//                                                  Element("Zn")});
//   const auto conf0 = cfg::GenerateFCC({8, 8, 8}, Element("Al"));
//   size_t Zn, Mg;
//   std::cout << "Mg Zn Energy" << std::endl;
// #pragma omp parallel for default(none) shared(conf0, a, std::cout) private(Zn, Mg)
//   for (Mg = 0; Mg <= 100; ++Mg) {
//     for (Zn = 0; Zn <= 100; ++Zn) {
//       if (Mg + Zn > 100) continue;
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

                                return 0;
                                }
