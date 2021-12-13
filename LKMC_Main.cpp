#include "Constants.hpp"
#include "VectorMatrix.hpp"
#include "Lattice.hpp"
#include "Atom.hpp"
#include "LatticeCluster.hpp"
#include "Config.h"
#include "KmcEvent.h"
#include "ChainKmcSimulation.h"

int main() {

  std::set<Element> ele_set{Element("Al"), Element("Mg"), Element("Zn")};
  auto conf = cfg::Config::ReadCfg("start.cfg");
  kmc::ChainKMCSimulation a(conf,
                            1e3,
                            1e5,
                            1e10,
                            ele_set,
                            0, 0, 0,
                            "kmc_parameters_bond.json");
  a.Simulate();

  // auto conf = cfg::Config::ReadCfg("forward.cfg");
  // pred::EnergyPredictorE0DECluster a("./kmc_parameters_cluster.json", conf, ele_set);
  // auto[Ea, dE] = a.GetBarrierAndDiffFromAtomIdPair(conf, {18, 23});
  // std::cout <<  Ea << ',' << dE << ',' << std::endl;

  // std::ifstream ifs("all_data_neb_results/barriers.txt", std::ifstream::in);
  // if (!ifs.is_open()) {
  //   std::cout << "Cannot open kmc_log.txt\n";
  //   return 1;
  // }
  // ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //
  // size_t vacancy_index, jump_index;
  // size_t ct = 0;
  // while (true) {
  //   if (ifs.eof() || ifs.bad()) {
  //     break;
  //   }
  //   ifs >> vacancy_index >> jump_index;
  //   ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //   std::string fn = "all_data_neb_results/config" + std::to_string(ct) + "/start.cfg";
  //   cfg::Config config = cfg::Config::ReadCfg(fn);
  //   pred::EnergyPredictorE0DEBond a("./kmc_parameters_bond.json", config, ele_set);
  //   auto[Ea, dE] = a.GetBarrierAndDiffFromAtomIdPair(config, {vacancy_index, jump_index});
  //   std::cout << ct << ',' << Ea << ',' << dE << ',' << std::endl;
  //   ct++;
  // }


  // auto map1 = pred::GetAverageClusterParametersMappingMMM(config);
  // auto map2 = pred::GetAverageClusterParametersMappingPeriodic(config);
  // // auto index_list = pred::GetSymmetricallySortedLatticeVectorMMM(config, {18, 23});
  // auto index_list = pred::GetSortedLatticeVectorPeriodic(config, 18);
  //
  // std::vector<Element> ele_vector{};
  // for (auto index: index_list) {
  //   ele_vector.push_back(config.GetElementAtLatticeID(index.GetId()));
  // }
  // auto res = pred::GetOneHotParametersFromMap(ele_vector, pred::GetOneHotEncodeHashmap(
  // ), 3, map2);

  // std::cout << '[';
  // for (auto i: res) {
  //   std::cout << i << ", ";
  // }
  // std::cout << ']';
  // cfg::Config config2 = cfg::Config::ReadCfg("backward.cfg");
  // auto map2 = pred::GetAverageClusterParametersMappingMMM(config2);
  // for (auto &i: map1) {
  //   for (auto &j: i) {
  //     std::sort(j.begin(), j.end());
  //   }
  //   std::sort(i.begin(), i.end());
  // }
  // for (const auto &i: map1) {
  //   for (const auto& j: i) {
  //     std::cout << '[';
  //     for (auto k: j) {
  //       std::cout << k << ' ';
  //     }
  //     std::cout << "]\n";
  //   }
  // }
  // std::cout << std::endl;

  // const auto &fal = config.GetFirstNeighborsAdjacencyList();
  // auto ll = config.GetLatticeVector();
  // for (size_t i = 0; i < fal.size(); ++i) {
  //   for (size_t j = 0; j < 12; ++j) {
  //     auto a = cfg::GetRelativeDistanceVectorLattice(ll[i], ll[fal[i][j]]);
  //     auto map = pred::GetAverageClusterParametersMappingMMM(config, {i, fal[i][j]});
  //
  //   }
  // }


  // size_t vacancy_site_index = cfg::GetVacancyLatticeIndex(config);
  // size_t neighbor_lattice_index = config.GetFirstNeighborsAdjacencyList()[vacancy_site_index][0];
  // std::cout << vacancy_site_index << '\t' << neighbor_lattice_index << '\n';
  // auto set = cfg::GetNeighborsLatticeIdSetOfJumpPair(
  //     config, {vacancy_site_index, neighbor_lattice_index});
  //
  // std::cout << set.size() << '\n';
  //
  // config.WriteCfg("same_as_start.cfg", false);
  // config.WriteCfg("same_as_0.cfg", true);

  return 0;
}
