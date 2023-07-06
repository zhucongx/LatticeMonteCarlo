/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/6/23 3:03 PM                                                            *
 **************************************************************************************************/

/*! \file  main.cpp
 *  \brief File for the main function.
 */

// #include "LatticeCluster.hpp"
// #include "AtomCluster.hpp"
// #include "ClusterExpansion.h"
// int main() {
//   std::cout << "Hello, World!" << std::endl;
//   Config cfg = Config::ReadCfg("test_small.cfg");
//   cfg.UpdateNeighborList({3.5, 4.8, 5.3, 5.9, 6.5, 7.1, 7.6});
//   auto tmp = InitializeLatticeSiteHashMap(cfg, 3, 3);
//   std::map<LatticeClusterType, size_t> tmp2(tmp.begin(), tmp.end());
//   for (auto [type, index] : tmp2) {
//     std::cout << type << " " << index << std::endl;
//   }
//   auto tmp3 = LatticeSiteHashMap(cfg, 3, 3);
//   std::map<LatticeClusterType, std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>>
//       tmp4(tmp3.begin(), tmp3.end());
//   for (auto [type, index] : tmp4) {
//     std::cout << type << " " << index.size() << std::endl;
//   }
//
//   std::cout << "Distance Order: " << cfg.GetDistanceOrder(0, 1) << std::endl;
//   // Config cfg = Config::ReadPoscar("POSCAR30.gz");
//   // cfg.SetPeriodicBoundaryCondition({false, false, false});
//   std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//   std::cout << "Neighbor list updated" << std::endl;
//   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//   std::cout << "Time difference = "
//             << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
//             << "[ms]" << std::endl;
//
//   std::vector<size_t> a(cfg.GetNumAtoms(), 0);
//   std::vector<size_t> b(cfg.GetNumAtoms(), 0);
//   std::vector<size_t> c(cfg.GetNumAtoms(), 0);
//   std::vector<size_t> d(cfg.GetNumAtoms(), 0);
//   std::vector<size_t> e(cfg.GetNumAtoms(), 0);
//   std::vector<size_t> f(cfg.GetNumAtoms(), 0);
//   std::vector<size_t> g(cfg.GetNumAtoms(), 0);
//
//   for (size_t i = 0; i < cfg.GetNumAtoms(); ++i) {
//     a[i] = cfg.GetNeighborLists()[0][i].size();
//     b[i] = cfg.GetNeighborLists()[1][i].size();
//     c[i] = cfg.GetNeighborLists()[2][i].size();
//     d[i] = cfg.GetNeighborLists()[3][i].size();
//     e[i] = cfg.GetNeighborLists()[4][i].size();
//     f[i] = cfg.GetNeighborLists()[5][i].size();
//     g[i] = cfg.GetNeighborLists()[6][i].size();
//   }
//   std::map<std::string, Config::VectorVariant>
//       auxiliary_lists{
//       {"num_nn1", a},
//       {"num_nn2", b},
//       {"num_nn3", c},
//       {"num_nn4", d},
//       {"num_nn5", e},
//       {"num_nn6", f},
//       {"num_nn7", g}
//   };
//   cfg.WriteXyzExtended("output.xyz.gz", auxiliary_lists, {});
//   cfg.WriteXyzExtended("output.xyz", auxiliary_lists, {});
//   cfg.WriteConfig("output.cfg.gz");
//   cfg.WriteConfig("output.cfg");
//   return 0;
// }

#include "Home.h"
int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  Parameter parameter(argc, argv);
  Print(parameter);
  Run(parameter);
  return 0;
}