#include "Traverse.h"

#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <omp.h>
using json = nlohmann::json;

namespace ansys {
static cfg::Config GetConfig(const std::string &config_type, size_t i) {
  cfg::Config config;
  if (config_type == "config") {
    config = cfg::Config::ReadConfig(std::to_string(i) + ".cfg.gz");
    config.ReassignLatticeVector();
  } else if (config_type == "map") {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", "map" + std::to_string(i) + ".txt");
  } else {
    throw std::invalid_argument("Unknown config type: " + config_type);
  }
  return config;
}

Traverse::Traverse(unsigned long long int initial_steps,
                   unsigned long long int increment_steps,
                   Element solvent_element,
                   std::set<Element> element_set,
                   size_t smallest_cluster_criteria,
                   size_t solvent_bond_criteria,
                   const std::string &predictor_filename,
                   std::string log_type,
                   std::string config_type)
    : initial_steps_(initial_steps),
      increment_steps_(increment_steps),
      final_number_(increment_steps),
      solvent_element_(solvent_element),
      element_set_(std::move(element_set)),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria),
      log_type_(std::move(log_type)),
      config_type_(std::move(config_type)),
      energy_predictor_(predictor_filename, element_set_),
      vacancy_migration_predictor_(predictor_filename, GetConfig(config_type_, 0), element_set_),
      energy_change_predictor_pair_site_(predictor_filename, GetConfig(config_type_, 0), element_set_) {
  std::string log_file_name;
  if (log_type_ == "kinetic_mc" || log_type_ == "kinetic_mc_old") {
    log_file_name = "kmc_log.txt";
  } else if (log_type_ == "canonical_mc") {
    log_file_name = "cmc_log.txt";
  } else {
    throw std::invalid_argument("Unknown log type: " + log_type_);
  }
  std::ifstream ifs(log_file_name, std::ifstream::in);
  if (!ifs.is_open()) {
    throw std::runtime_error("Cannot open " + log_file_name);
  }
  unsigned long long step_number;
  double energy{}, time{}, temperature{};
  while (ifs.peek() != '0') {
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  while (true) {
    if (ifs.eof() || ifs.bad()) {
      break;
    }
    if (ifs.peek() == '#') {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    if (log_type_ == "kinetic_mc") {
      ifs >> step_number >> time >> temperature >> energy;
    } else if (log_type_ == "kinetic_mc_old") {
      ifs >> step_number >> time >> energy;
    } else if (log_type_ == "canonical_mc") {
      ifs >> step_number >> temperature >> energy;
    } else {
      throw std::invalid_argument("Unknown log type: " + log_type_);
    }

    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (step_number >= initial_steps_ && (step_number - initial_steps_) % increment_steps == 0) {
      filename_energy_hashset_[step_number] = energy;
      filename_temperature_hashset_[step_number] = temperature;
      filename_time_hashset_[step_number] = time;
      final_number_ = step_number;
    }
  }
#pragma omp parallel default(none) shared(std::cout)
  {
#pragma omp master
    { std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl; }
  }
}

Traverse::~Traverse() = default;

void Traverse::RunAnsys() const {
  const auto chemical_potential = energy_predictor_.GetChemicalPotential(solvent_element_);
  json ansys_info_array = json::array();
#pragma omp parallel for default(none) schedule(static, 1) shared(ansys_info_array, chemical_potential, std::cout)
  for (unsigned long long i = initial_steps_; i <= final_number_; i += increment_steps_) {
#pragma omp critical
    {
      std::cout << i << " / " << final_number_ << " " << std::fixed << std::setprecision(2)
                << static_cast<double>(i) / static_cast<double>(final_number_) * 100 << "%" << std::endl;
    }
    cfg::Config config = GetConfig(config_type_, i);

    // basic information
    json ansys_info = json::object();
    std::map<std::string, cfg::Config::ValueVariant> global_list;
    ansys_info["steps"] = i;
    global_list["steps"] = i;

    ansys_info["time"] = filename_time_hashset_.at(i);
    global_list["time"] = filename_time_hashset_.at(i);

    ansys_info["temperature"] = filename_temperature_hashset_.at(i);
    global_list["temperature"] = filename_temperature_hashset_.at(i);

    ansys_info["energy"] = filename_energy_hashset_.at(i);
    global_list["energy"] = filename_energy_hashset_.at(i);

    // cluster information
    auto [cluster_json, auxiliary_lists] = SoluteCluster(config,
                                                         solvent_element_,
                                                         element_set_,
                                                         smallest_cluster_criteria_,
                                                         solvent_bond_criteria_,
                                                         energy_predictor_,
                                                         chemical_potential)
                                               .GetClustersInfo();
    ansys_info["clusters"] = cluster_json;

    // sro information
    ShortRangeOrder short_range_order(config, element_set_);
    const auto sro1 = short_range_order.FindWarrenCowley(1);
    ansys_info["warren_cowley"]["first"] = sro1;
    for (const auto &[pair, value]: sro1) {
      global_list["sro1_" + pair] = value;
    }

    const auto sro2 = short_range_order.FindWarrenCowley(2);
    ansys_info["warren_cowley"]["second"] = sro2;
    for (const auto &[pair, value]: sro2) {
      global_list["sro2_" + pair] = value;
    }

    const auto sro3 = short_range_order.FindWarrenCowley(3);
    ansys_info["warren_cowley"]["third"] = sro3;
    for (const auto &[pair, value]: sro3) {
      global_list["sro3_" + pair] = value;
    }

    // vacancy information
    auto vacancy_lattice_id = config.GetVacancyLatticeId();
    static auto convert = [](const std::set<Element> &element_set_base,
                             const std::map<Element, size_t> &element_map) -> std::map<std::string, size_t> {
      std::map<std::string, size_t> string_map;
      for (const auto &element: element_set_base) {
        string_map[element.GetString()] = 0;
      }
      for (const auto &[element, value]: element_map) {
        string_map[element.GetString()] = value;
      }
      return string_map;
    };
    ansys_info["vac_local"]["first"] = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 1));
    ansys_info["vac_local"]["second"] = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 2));
    ansys_info["vac_local"]["third"] = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 3));

    // binding energy
    const auto exit_time = ExitTime(config,
                                    solvent_element_,
                                    filename_temperature_hashset_.at(i),
                                    vacancy_migration_predictor_,
                                    energy_change_predictor_pair_site_,
                                    chemical_potential);
    const auto binding_energy = exit_time.GetBindingEnergy();
    auxiliary_lists["binding_energy"] = binding_energy;

    ansys_info["vac_local_binding_energy"] = binding_energy[config.GetVacancyLatticeId()];
    global_list["vac_local_binding_energy"] = binding_energy[config.GetVacancyLatticeId()];

    for (auto &cluster_info: ansys_info["clusters"]) {
      std::vector<double> binding_energy_list;
      for (const auto &atom_id: cluster_info["cluster_atom_id_list"]) {
        binding_energy_list.push_back(binding_energy[atom_id]);
      }
      cluster_info.erase("cluster_atom_id_list");
      cluster_info["vacancy_binding_energy"] =
          *std::min_element(binding_energy_list.begin(), binding_energy_list.end());
    }

    // const auto profile_energy = exit_time.GetProfileEnergy();
    // auxiliary_lists["profile_energy"] = profile_energy;

    // exit time
    // auto [barrier_lists, exit_times] =exit_time.GetBarrierListAndExitTime();
    // auxiliary_lists["barrier_lists"] = barrier_lists;
    // auxiliary_lists["exit_times"] = exit_times;

    boost::filesystem::create_directories("ansys");
    config.WriteExtendedXyz("ansys/" + std::to_string(i) + ".xyz.gz", auxiliary_lists, global_list);
#pragma omp critical
    { ansys_info_array.push_back(ansys_info); }
  }
  std::cout << "Finished. Sorting..." << std::endl;
  std::sort(ansys_info_array.begin(), ansys_info_array.end(), [](const json &lhs, const json &rhs) {
    return lhs["steps"] < rhs["steps"];
  });
  std::ofstream ofs("ansys_info.json.gz", std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  fos.push(boost::iostreams::gzip_compressor());
  fos.push(ofs);
  fos.precision(16);
  fos << ansys_info_array.dump(2) << std::endl;
  std::cout << "Done..." << std::endl;
}

void Traverse::RunReformat() const {
  for (unsigned long long i = 0; i <= final_number_; i += increment_steps_) {
    std::cout << i << " / " << final_number_ << std::endl;
    if (config_type_ == "map") {
      auto config = cfg::Config::ReadMap("lattice.txt", "element.txt", "map" + std::to_string(i) + ".txt");
      config.WriteConfig(std::to_string(i) + ".cfg.gz");
    } else if (config_type_ == "config") {
      auto config = cfg::Config::ReadConfig(std::to_string(i) + ".cfg");
      config.ReassignLatticeVector();
      if (i == 0) {
        config.WriteLattice("lattice.txt");
        config.WriteElement("element.txt");
      }
      config.WriteMap("map" + std::to_string(i) + ".txt");
    } else {
      throw std::invalid_argument("Unknown config type: " + config_type_);
    }
  }
}
}    // namespace ansys