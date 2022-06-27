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

  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << std::setprecision(16)
  //           << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();


  const pred::TotalEnergyPredictor energy_estimator("quartic_coefficients.json",
                                                    {Element("Al"),
                                                     Element("Mg"),
                                                     Element("Zn"),
                                                    });

  std::cout << "Energy 444: " << energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("444_c.cfg")) - energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("444_f.cfg")) << std::endl;
  std::cout << "Energy 555: " << energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("555_c.cfg")) - energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("555_f.cfg")) << std::endl;

  pred::StateChangePredictor a("quartic_coefficients.json",
                               cfg::Config::ReadCfg("444_c.cfg"),
                               {Element("Al"),
                                Element("Mg"),
                                Element("Zn")});
  std::cout << "Energy change 444: "
            << a.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("444_c.cfg"), {1, 168})
            << std::endl;
  std::cout << "Energy change 444: "
            << a.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("444_c.cfg"), {168, 1})
            << std::endl;

  pred::StateChangePredictor b("quartic_coefficients.json",
                               cfg::Config::ReadCfg("555_c.cfg"),
                               {Element("Al"),
                                Element("Mg"),
                                Element("Zn")});
  std::cout << "Energy change 555: "
            << b.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("555_c.cfg"), {1, 349})
            << std::endl;
  std::cout << "Energy change 555: "
            << b.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("555_c.cfg"), {349, 1})
            << std::endl;

  auto conf444 = cfg::Config::ReadCfg("444_b.cfg");
  const pred::VacancyMigrationPredictorQuartic
      energy_predictor444("quartic_coefficients.json", conf444,
                          {Element("Al"),
                           Element("Mg"),
                           Element("Zn")});
  auto pair444 = energy_predictor444.GetBarrierAndDiffFromAtomIdPair(conf444, {0, 61});
  std::cout << "Energy barrier 444: " << pair444.first << " Energy change 444: " << pair444.second
            << std::endl;

  auto conf555 = cfg::Config::ReadCfg("555_b.cfg");
  const pred::VacancyMigrationPredictorQuartic
      energy_predictor555("quartic_coefficients.json", conf555,
                          {Element("Al"),
                           Element("Mg"),
                           Element("Zn")});
  auto pair555 = energy_predictor555.GetBarrierAndDiffFromAtomIdPair(conf555, {0, 97});
  std::cout << "Energy barrier 555: " << pair555.first << " Energy change 555: " << pair555.second
            << std::endl;

  pred::StateChangePredictor c("quartic_coefficients.json",
                               cfg::Config::ReadCfg("444_b.cfg"),
                               {Element("Al"),
                                Element("Mg"),
                                Element("Zn")});
  std::cout << "Energy change 444: "
            << c.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("444_b.cfg"), {61, 0})
            << std::endl;
  std::cout << "Energy change 444: "
            << c.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("444_b.cfg"), {0, 61})
            << std::endl;

  pred::StateChangePredictor d("quartic_coefficients.json",
                               cfg::Config::ReadCfg("555_b.cfg"),
                               {Element("Al"),
                                Element("Mg"),
                                Element("Zn")});
  std::cout << "Energy change 555: "
            << d.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("555_b.cfg"), {0, 97})
            << std::endl;
  std::cout << "Energy change 555: "
            << d.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("555_b.cfg"), {97, 0})
            << std::endl;

  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < 100; ++i)
    auto aa = d.GetDiffFromAtomIdPair(cfg::Config::ReadCfg("555_b.cfg"), {97, 0});
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

//   pred::TotalEnergyPredictor a("quartic_coefficients.json",
//                                std::set<Element>{Element("Al"), Element("Mg"),
//                                                  Element("Zn")});
//   const auto conf0 = cfg::GenerateFCC(
//       4.046, {8, 8, 8}, Element("Al"));
//   size_t Zn, Mg;
//   std::cout << "Mg Zn Energy" << std::endl;
// #pragma omp parallel for default(none) shared(conf0, a, std::cout) private(Zn, Mg)
//   for (Mg = 0; Mg <= 100; ++Mg) {
//     for (Zn = 0; Zn <= 100; ++Zn) {
//       if (Mg + Zn > 100) continue;
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
