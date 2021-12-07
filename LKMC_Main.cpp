#include "Constants.hpp"
#include "VectorMatrix.hpp"
#include "Lattice.hpp"
#include "Atom.hpp"
#include "LatticeCluster.hpp"
#include "Config.h"
#include "ClusterExpansionPredictor.h"
bool operator<(const std::vector<int> &lhs, const std::vector<int> &rhs) {
  for (size_t i =0 ; i <lhs.size(); ++i) {
    if (lhs[i] < rhs[i]) return true;
    return false;
  }
}

int main() {
  class ClusterExpansionPredictor;
  cfg::Config config = cfg::Config::ReadCfg("forward.cfg");
  auto map1 = pred::GetAverageClusterParametersMappingPeriodic(config, 18);

  cfg::Config config2 = cfg::Config::ReadCfg("backward.cfg");
  auto map2 = pred::GetAverageClusterParametersMappingPeriodic(config2, 18);
  for (auto &i: map1) {
    for (auto &j: i) {
      std::sort(j.begin(), j.end());
    }
    std::sort(i.begin(), i.end());
  }
  for (const auto &i: map1) {
    for (const auto& j: i) {
      std::cout << '[';
      for (auto k: j) {
        std::cout << k << ' ';
      }
      std::cout << "]\n";
    }
  }
  std::cout << std::endl;

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
  // auto set = cfg::GetFirstAndSecondThirdNeighborsLatticeIdSetOfJumpPair(
  //     config, {vacancy_site_index, neighbor_lattice_index});
  //
  // std::cout << set.size() << '\n';
  //
  // config.WriteCfg("same_as_start.cfg", false);
  // config.WriteCfg("same_as_0.cfg", true);

  return 0;
}
