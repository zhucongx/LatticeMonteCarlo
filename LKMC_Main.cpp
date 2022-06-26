#include <chrono>
#include "Home.h"
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
  // }

  // ansys::Iterator test(0, 1e4,
  //                      Element("Al"),
  //                      {Element("Al"),
  //                       Element("Mg"),
  //                       Element("Zn")},
  //                      4, "quartic_coefficients.json");
  // // test.SerialRunReformat();
  // test.SerialRunCluster();


  // auto t1 = std::chrono::high_resolution_clock::now();
  // const pred::TotalEnergyEstimator energy_estimator("quartic_coefficients.json",
  //                                                   {Element("Al"),
  //                                                    Element("Mg"),
  //                                                    Element("Zn")});
  // std::cout << "Energy 333: " << energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("333_f.cfg")) - energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("333_c.cfg")) << std::endl;
  // std::cout << "Energy 444: " << energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("444_f.cfg")) - energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("444_c.cfg")) << std::endl;
  // std::cout << "Energy 555: " << energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("555_f.cfg")) - energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("555_c.cfg")) << std::endl;
  // std::cout << "Energy 303030: " << energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("303030_f.cfg")) - energy_estimator.GetEnergy(
  //     cfg::Config::ReadCfg("303030_c.cfg")) << std::endl;
  // auto conf333 = cfg::Config::ReadCfg("333_b.cfg");
  // const pred::VacancyMigrationPredictorQuartic
  //     energy_predictor333("quartic_coefficients.json", conf333,
  //                         {Element("Al"),
  //                          Element("Mg"),
  //                          Element("Zn")});
  // auto pair333 = energy_predictor333.GetBarrierAndDiffFromAtomIdPair(conf333, {0, 33});
  // std::cout << "Energy barrier 333: " << pair333.first << " Energy change 333: " << pair333.second
  //           << std::endl;
  //
  // auto conf444 = cfg::Config::ReadCfg("444_b.cfg");
  // const pred::VacancyMigrationPredictorQuartic
  //     energy_predictor444("quartic_coefficients.json", conf444,
  //                         {Element("Al"),
  //                          Element("Mg"),
  //                          Element("Zn")});
  // auto pair444 = energy_predictor444.GetBarrierAndDiffFromAtomIdPair(conf444, {0, 61});
  // std::cout << "Energy barrier 444: " << pair444.first << " Energy change 444: " << pair444.second
  //           << std::endl;
  //
  // auto conf555 = cfg::Config::ReadCfg("555_b.cfg");
  // const pred::VacancyMigrationPredictorQuartic
  //     energy_predictor555("quartic_coefficients.json", conf555,
  //                         {Element("Al"),
  //                          Element("Mg"),
  //                          Element("Zn")});
  // auto pair555 = energy_predictor555.GetBarrierAndDiffFromAtomIdPair(conf555, {0, 97});
  // std::cout << "Energy barrier 555: " << pair555.first << " Energy change 555: " << pair555.second
  //           << std::endl;
  //
  // auto conf303030 = cfg::Config::ReadCfg("303030_b.cfg");
  // const pred::VacancyMigrationPredictorQuartic
  //     energy_predictor303030("quartic_coefficients.json", conf303030,
  //                            {Element("Al"),
  //                             Element("Mg"),
  //                             Element("Zn")});
  // auto pair303030 = energy_predictor303030.GetBarrierAndDiffFromAtomIdPair(conf303030, {0, 3597});
  // std::cout << "Energy barrier 303030: " << pair303030.first << " Energy change 303030: "
  //           << pair303030.second
  //           << std::endl;
  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << std::setprecision(16)
  //           << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();


//   pred::TotalEnergyEstimator a("quartic_coefficients.json",
//                                std::set<Element>{Element("Al"), Element("Mg"),
//                                                  Element("Zn")});
//   const auto conf0 = cfg::GenerateFCC(
//       4.046, {10, 10, 10}, Element("Al"));
//   size_t Zn, Mg;
//   std::cout << "Mg Zn Energy" << std::endl;
// #pragma omp parallel for default(none) shared(conf0, a, std::cout) private(Zn, Mg)
//   for (Mg = 0; Mg <= 200; ++Mg) {
//     for (Zn = 0; Zn <= 200; ++Zn) {
//       if (Mg + Zn > 200) continue;
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
