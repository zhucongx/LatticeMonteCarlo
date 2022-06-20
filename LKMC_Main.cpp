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
  //                       4, 3, "quartic_coefficients.json");
  // test.SerialRunCluster();
  // auto conf0 = cfg::GenerateFCC(4.046, {30, 30, 30}, Element("Al"));
  // conf0.WriteCfg("303030.cfg", true);
  auto t1 = std::chrono::high_resolution_clock::now();
  const pred::EnergyEstimator energy_estimator("quartic_coefficients.json",
                                               {Element("Al"),
                                                Element("Mg"),
                                                Element("Zn")});
  std::cout << "Energy 333: " << energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("333_f.cfg")) - energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("333_c.cfg")) << std::endl;

  std::cout << "Energy 444: " << energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("444_f.cfg")) - energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("444_c.cfg")) << std::endl;

  std::cout << "Energy 555: " << energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("555_f.cfg")) - energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("555_c.cfg")) << std::endl;
  std::cout << "Energy 303030: " << energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("303030_f.cfg")) - energy_estimator.GetEnergy(
      cfg::Config::ReadCfg("303030_c.cfg")) << std::endl;

  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << std::endl;
  // auto cluster_finder = ansys::ClustersFinder("22100000.cfg", Element("Al"),
  //                                             4, 3, energy_estimator);
  // auto result = cluster_finder.FindClustersAndOutput();
  // for (auto [cluster, energy]: result) {
  //   std::cout << energy << '\n';
  // }

  // auto conf0 = cfg::GenerateFCC(4.046, {4,4,4}, Element("Al"));
  //
  // conf0.WriteCfg("444.cfg", true);
  // std::cout << conf0.GetStateHash() << std::endl;
  //
  // auto conf1 = cfg::Config::ReadCfg("start1.cfg");
  // conf1.WriteCfg("start2.cfg", true);
  // std::cout << conf1.GetStateHash() << std::endl;
  //
  // conf1.WriteLattice("start_lattice.txt");
  // conf1.WriteElement("start_element.txt");
  // conf1.WriteMap("start_map.txt");
  //
  // auto conf2 = cfg::Config::ReadMap("start_lattice.txt", "start_element.txt", "start_map.txt");
  // conf2.WriteCfg("start3.cfg", true);
  // std::cout << conf2.GetStateHash() << std::endl;

  // std::set<Element> ele_set{Element("Al"), Element("Mg"),
  //                           Element("Zn")};
  // auto conf = cfg::Config::ReadCfg("start.cfg");
  // pred::EnergyPredictorQuarticLru a("./kmc_parameters_quartic.json",
  //                                 conf, ele_set, 100000);
  // auto t1 = std::chrono::high_resolution_clock::now();
  // auto [Ea, dE] = a.GetBarrierAndDiffFromAtomIdPair(conf, {82, 83});
  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << std::setprecision(16)
  //           << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << ", "
  //           << Ea << ", " << dE << std::endl;
  // for (int i = 0; i < 1000; ++i) {
  //   auto pair = a.GetBarrierAndDiffFromAtomIdPair(conf, {82, 83});
  //   Ea = pair.first;
  //   dE = pair.second;
  // }
  // auto t3 = std::chrono::high_resolution_clock::now();
  // std::cout << std::setprecision(16)
  //           << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << ", "
  //           << Ea << ", " << dE << std::endl;


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
  // std::set<Element> ele_set1{Element("Al"),
  //                           Element("Mg"),
  //                           Element("Zn"),
  //                           Element("X"),
  //                           Element("pZn"),
  //
  // };
  //
  // auto conf1 = cfg::Config::ReadCfg("forward.cfg");
  //
  // auto mapping_state_ = pred::GetClusterParametersMappingState(conf1);
  // auto ea = pred::GetEncodesFromMapState(conf1, {18, 23}, pred::InitializeClusterHashMap(ele_set1),
  //                                        mapping_state_);
  //
  //

  // auto one_hot_encode_hashmap = pred::InitializeClusterHashMap(ele_set);
  //
  // auto s_sorted_lattice_vector_periodic = pred::GetSortedLatticeVectorPeriodic(conf, 18);
  // std::vector<size_t> s_lattice_id_vector;
  // std::transform(s_sorted_lattice_vector_periodic.begin(), s_sorted_lattice_vector_periodic.end(),
  //                std::back_inserter(s_lattice_id_vector),
  //                [](const auto &lattice) { return lattice.GetId(); });
  // auto e_sorted_lattice_vector_periodic = pred::GetSortedLatticeVectorPeriodic(conf, 23);
  // std::vector<size_t> e_lattice_id_vector;
  // std::transform(e_sorted_lattice_vector_periodic.begin(), e_sorted_lattice_vector_periodic.end(),
  //                std::back_inserter(e_lattice_id_vector),
  //                [](const auto &lattice) { return lattice.GetId(); });
  //
  // std::vector<Element> element_vector_start{}, element_vector_end{};
  // for (auto index: s_lattice_id_vector) {
  //   auto this_element = conf.GetElementAtLatticeId(index);
  //   if (this_element == ElementName::X) {
  //     element_vector_start.push_back(conf.GetElementAtLatticeId(23));
  //     continue;
  //   }
  //   element_vector_start.push_back(this_element);
  // }
  // auto start_encode = pred::GetClusterEncode(element_vector_start,
  //                                      mapping_periodic_,
  //                                            one_hot_encode_hashmap);
  //
  //
  // for (auto index: e_lattice_id_vector) {
  //   auto this_element = conf.GetElementAtLatticeId(index);
  //   if (this_element == ElementName::X) {
  //     element_vector_end.push_back(conf.GetElementAtLatticeId(23));
  //     continue;
  //   }
  //   element_vector_end.push_back(this_element);
  // }
  // auto end_encode = GetClusterEncode(element_vector_end,
  //                                    mapping_periodic_,
  //                                    one_hot_encode_hashmap);
  //
  // for (size_t i = 0; i < end_encode.size(); ++i) {
  //   std::cout << end_encode[i] - start_encode[i] << ", ";
  // }
  // std::cout << "\n";
  //
  // auto index_list = pred::GetSymmetricallySortedLatticeVectorMMM(config, {18, 23});
  // auto index_list = pred::GetSortedLatticeVectorPeriodic(config, 18);
  //
  // std::vector<Element> ele_vector{};
  // for (auto index: index_list) {
  //   ele_vector.push_back(config.GetElementAtLatticeID(index.GetId()));
  // }
  //
  // auto res = pred::GetOneHotParametersFromMapPeriodic(conf, ele_set, map2);
  //
  // // std::cout << '[';
  // for (auto i: res) {
  //   std::cout << i << ", ";
  // }
  // std::cout << ']';
  // cfg::Config config2 = cfg::Config::ReadCfg("backward.cfg");
  // auto map2 = pred::GetAverageClusterParametersMappingMMM(config2);
  // for (auto &i: map2) {
  //   for (auto &j: i) {
  //     std::sort(j.begin(), j.end());
  //   }
  //   std::sort(i.begin(), i.end());
  // }
  // for (const auto &i: map2) {
  //   for (const auto& j: i) {
  //     for (auto k: j) {
  //       std::cout << k << ',';
  //     }
  //     std::cout << "\n";
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
