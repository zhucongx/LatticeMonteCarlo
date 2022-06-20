#include "EnergyEstimator.h"
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
EnergyEstimator::EnergyEstimator(const std::string &predictor_filename,
                                 std::set<Element> type_set)
    : type_set_(std::move(type_set)) {
  auto type_set_copy(type_set_);
  type_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(type_set_copy);

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());

  // for (const auto &cluster_count: ordered) {
  //   sorted_cluster_type_vector.push_back( cluster_count.first);
  // }

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
}
EnergyEstimator::~EnergyEstimator() = default;
std::vector<double> EnergyEstimator::GetEncodeFast(const cfg::Config &config) const {
  auto cluster_hashmap(initialized_cluster_hashmap_);
  const std::vector<size_t> cluster_counter{256,
                                            3072,
                                            1536,
                                            6144,
                                            12288,
                                            6144,
                                            12288,
                                            6144,
                                            12288,
                                            12288,
                                            12288};
  for (size_t atom_id1 = 0; atom_id1 < config.GetNumAtoms(); ++atom_id1) {
    const size_t lattice_id1 = config.GetLatticeIdFromAtomId(atom_id1);
    Element element1 = config.GetAtomVector()[atom_id1].GetElement();
    cluster_hashmap[cfg::ElementCluster(0, element1)]++;
    for (size_t lattice_id2: config.GetFirstNeighborsAdjacencyList()[lattice_id1]) {
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(1, element1, element2)]++;
      for (size_t lattice_id3: config.GetFirstNeighborsAdjacencyList()[lattice_id2]) {
        Element element3 = config.GetElementAtLatticeId(lattice_id3);
        if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice_id1].begin(),
                      config.GetFirstNeighborsAdjacencyList()[lattice_id1].end(),
                      lattice_id3)
            != config.GetFirstNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(4, element1, element2, element3)]++;
        } else if (std::find(config.GetSecondNeighborsAdjacencyList()[lattice_id1].begin(),
                             config.GetSecondNeighborsAdjacencyList()[lattice_id1].end(),
                             lattice_id3)
            != config.GetSecondNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(5, element1, element2, element3)]++;
        } else if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
                             config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
                             lattice_id3)
            != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(6, element1, element2, element3)]++;
        }
      }
      for (size_t lattice_id3: config.GetSecondNeighborsAdjacencyList()[lattice_id2]) {
        Element element3 = config.GetElementAtLatticeId(lattice_id3);
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
                      lattice_id3)
            != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(7, element1, element2, element3)]++;
        }
      }
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(8, element1, element2, element3)]++;
      //   }
      // }
    }
    for (size_t lattice_id2: config.GetSecondNeighborsAdjacencyList()[lattice_id1]) {
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(2, element1, element2)]++;
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(9, element1, element2, element3)]++;
      //   }
      // }
    }
    for (size_t lattice_id2: config.GetThirdNeighborsAdjacencyList()[lattice_id1]) {
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(3, element1, element2)]++;
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(10, element1, element2, element3)]++;
      //   }
      // }
    }
  }

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> energy_encode;
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto count_bond = static_cast<double>(cluster_hashmap.at(cluster));
    auto total_bond = static_cast<double>(cluster_counter[static_cast<size_t>(cluster.GetLabel())]);
    energy_encode.push_back(count_bond / total_bond);
  }
  return energy_encode;
}
std::vector<double> EnergyEstimator::GetEncodeFastOmp(const cfg::Config &config) const {
  auto cluster_hashmap(initialized_cluster_hashmap_);
  const std::vector<size_t> cluster_counter{256,
                                            3072,
                                            1536,
                                            6144,
                                            12288,
                                            6144,
                                            12288,
                                            6144,
                                            12288,
                                            12288,
                                            12288};
#pragma omp parallel default(none) shared(config, cluster_hashmap, cluster_counter)
  {
#pragma omp for
    for (size_t atom_id1 = 0; atom_id1 < config.GetNumAtoms(); ++atom_id1) {
      const size_t lattice_id1 = config.GetLatticeIdFromAtomId(atom_id1);
      Element element1 = config.GetAtomVector()[atom_id1].GetElement();
      auto count0 =
          cluster_hashmap.find(cfg::ElementCluster(0, element1));
#pragma omp critical
      {
        count0->second++;
      }
      for (size_t lattice_id2: config.GetFirstNeighborsAdjacencyList()[lattice_id1]) {
        Element element2 = config.GetElementAtLatticeId(lattice_id2);
        auto count1 =
            cluster_hashmap.find(cfg::ElementCluster(1, element1, element2));
#pragma omp critical
        {
          count1->second++;
        }
        for (size_t lattice_id3: config.GetFirstNeighborsAdjacencyList()[lattice_id2]) {
          Element element3 = config.GetElementAtLatticeId(lattice_id3);
          if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice_id1].begin(),
                        config.GetFirstNeighborsAdjacencyList()[lattice_id1].end(),
                        lattice_id3)
              != config.GetFirstNeighborsAdjacencyList()[lattice_id1].end()) {
            auto count4 =
                cluster_hashmap.find(cfg::ElementCluster(4, element1, element2, element3));
#pragma omp critical
            {
              count4->second++;
            }
          } else if (std::find(config.GetSecondNeighborsAdjacencyList()[lattice_id1].begin(),
                               config.GetSecondNeighborsAdjacencyList()[lattice_id1].end(),
                               lattice_id3)
              != config.GetSecondNeighborsAdjacencyList()[lattice_id1].end()) {
            auto count5 =
                cluster_hashmap.find(cfg::ElementCluster(5, element1, element2, element3));
#pragma omp critical
            {
              count5->second++;
            }
          } else if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
                               config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
                               lattice_id3)
              != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
            auto count6 =
                cluster_hashmap.find(cfg::ElementCluster(6, element1, element2, element3));
#pragma omp critical
            {
              count6->second++;
            }
          }
        }
        for (size_t lattice_id3: config.GetSecondNeighborsAdjacencyList()[lattice_id2]) {
          Element element3 = config.GetElementAtLatticeId(lattice_id3);
          if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
                        config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
                        lattice_id3)
              != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
            auto count7 =
                cluster_hashmap.find(cfg::ElementCluster(7, element1, element2, element3));
#pragma omp critical
            {
              count7->second++;
            }
          }
        }
//         for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
//           Element element3 = config.GetElementAtLatticeId(lattice_id3);
//           if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
//                         config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
//                         lattice_id3)
//               != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
//             auto count8 =
//                 cluster_hashmap.find(cfg::ElementCluster(8, element1, element2, element3));
// #pragma omp critical
//             {
//               count8->second++;
//             }
//           }
//         }
      }
      for (size_t lattice_id2: config.GetSecondNeighborsAdjacencyList()[lattice_id1]) {
        Element element2 = config.GetElementAtLatticeId(lattice_id2);
        auto count2 =
            cluster_hashmap.find(cfg::ElementCluster(2, element1, element2));
#pragma omp critical
        {
          count2->second++;
        }
//         for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
//           Element element3 = config.GetElementAtLatticeId(lattice_id3);
//           if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
//                         config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
//                         lattice_id3)
//               != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
//             auto count9 =
//                 cluster_hashmap.find(cfg::ElementCluster(9, element1, element2, element3));
// #pragma omp critical
//             {
//               count9->second++;
//             }
//           }
//         }
      }
      for (size_t lattice_id2: config.GetThirdNeighborsAdjacencyList()[lattice_id1]) {
        Element element2 = config.GetElementAtLatticeId(lattice_id2);
        auto count3 =
            cluster_hashmap.find(cfg::ElementCluster(3, element1, element2));
#pragma omp critical
        {
          count3->second++;
        }
//         for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
//           Element element3 = config.GetElementAtLatticeId(lattice_id3);
//           if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
//                         config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
//                         lattice_id3)
//               != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
//             auto count10 =
//                 cluster_hashmap.find(cfg::ElementCluster(10, element1, element2, element3));
// #pragma omp critical
//             {
//               count10->second++;
//             }
//           }
//         }
      }
    }
  }
  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> energy_encode;
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto count_bond = static_cast<double>(cluster_hashmap.at(cluster));
    auto total_bond = static_cast<double>(cluster_counter[static_cast<size_t>(cluster.GetLabel())]);
    energy_encode.push_back(count_bond / total_bond);
  }
  return energy_encode;
}
double EnergyEstimator::GetEnergy(const cfg::Config &config) const {
  auto encode = GetEncodeFast(config);
  double E = 0;
  const size_t cluster_size = base_theta_.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    E += base_theta_[i] * encode[i];
  }
  return E;
}

} // namespace pred