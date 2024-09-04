#include "Home.h"

int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
  return 0;
}
// int main(int argc, char *argv[]) {
//   std::set<Element> element_set{Element("Al"), Element("Mg"), Element("Zn"), Element("X")};
//   const auto energy_estimator = pred::EnergyPredictor("quartic_coefficients.json", element_set);
//
//   auto cfg = cfg::Config::ReadConfig("./gp_cubic.cfg");
//   std::cout << energy_estimator.GetEnergy(cfg) << std::endl;
//
//   cfg = cfg::Config::ReadConfig("./lowest_energy.cfg");
//   std::cout << energy_estimator.GetEnergy(cfg) << std::endl;
//
//   cfg = cfg::Config::ReadConfig("./POSCAR1.cfg");
//   std::cout <<std::setprecision (15) << energy_estimator.GetEnergy(cfg) << std::endl;
//
//   cfg = cfg::Config::ReadConfig("./POSCAR2.cfg");
//   std::cout << std::setprecision (15) << energy_estimator.GetEnergy(cfg) << std::endl;
//
//   cfg = cfg::Config::ReadConfig("./POSCAR3.cfg");
//   std::cout << std::setprecision (15) << energy_estimator.GetEnergy(cfg) << std::endl;
// }

// int i = 800;
// while (i >= 275) {
//   auto cfg = cfg::Config::ReadConfig("/Users/zhucongx/Research/GOALI/mc_new/cmc_new/fetch_small/end_" + std::to_string(i) +
//                                      "K.cfg.gz");
//   cfg.ReassignLatticeVector();
//   // std::cout << i << "\t";
//   std::cout << std::setprecision(16)<< energy_estimator.GetEnergy(cfg) << ", " << std::endl;
//   i -= 25;
// }
// auto solvent_config = cfg::GenerateFCC({10, 10, 10}, Element("Al"));
// solvent_config.SetAtomElementTypeAtLattice(0, Element("X"));
// auto energy_predictor = pred::EnergyChangePredictorPairSite(
//     "quartic_coefficients_AlMgSnZn_2.json", solvent_config, {Element("Al"), Element("Mg"), Element("Sn"), Element("Zn")});
//
// auto neighbor_lattice1 = solvent_config.GetFirstNeighborsAdjacencyList()[0][0];
// auto neighbor_lattice2 = solvent_config.GetSecondNeighborsAdjacencyList()[0][0];
// auto neighbor_lattice3 = solvent_config.GetThirdNeighborsAdjacencyList()[0][0];
//
//
// solvent_config.SetAtomElementTypeAtLattice(2000, Element("Sn"));
// std::cout << "energy_change1 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice1})
//           << std::endl;
// std::cout << "energy_change2 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice2})
//           << std::endl;
// std::cout << "energy_change3 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice3})
//           << std::endl;
//
// solvent_config.SetAtomElementTypeAtLattice(2000, Element("Zn"));
// std::cout << "energy_change1 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice1})
//           << std::endl;
// std::cout << "energy_change2 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice2})
//           << std::endl;
// std::cout << "energy_change3 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice3})
//           << std::endl;
//
// solvent_config.SetAtomElementTypeAtLattice(2000, Element("Mg"));
// std::cout << "energy_change1 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice1})
//           << std::endl;
// std::cout << "energy_change2 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice2})
//           << std::endl;
// std::cout << "energy_change3 = " << energy_predictor.GetDeFromLatticeIdPair(solvent_config, {2000, neighbor_lattice3})
//           << std::endl;
// }