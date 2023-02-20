#include "Iterator.h"

#include <utility>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <experimental/iterator>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace ansys {

Iterator::Iterator(unsigned long long int initial_steps,
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
      energy_estimator_(predictor_filename, element_set_),
      log_type_(std::move(log_type)),
      config_type_(std::move(config_type)) {
  std::string log_file_name;
  if (log_type_ == "kinetic_mc") {
    log_file_name = "kmc_log.txt";
  } else if (log_type_ == "canonical_mc") {
    log_file_name = "cmc_log.txt";
  } else {
    std::cerr << "Unknown log type: " << log_type_ << std::endl;
  }
  std::ifstream ifs(log_file_name, std::ifstream::in);
  if (!ifs.is_open()) {
    std::cerr << "Cannot open " << log_file_name << "\n";
    return;
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
    } else if (log_type_ == "canonical_mc") {
      ifs >> step_number >> temperature >> energy;
    } else {
      std::cerr << "Unknown log type: " << log_type_ << std::endl;
    }

    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (step_number >= initial_steps_ && (step_number - initial_steps_) % increment_steps == 0) {
      filename_energy_hashset_[step_number] = energy;
      filename_temperature_hashset_[step_number] = temperature;
      filename_time_hashset_[step_number] = time;

      final_number_ = step_number;
    }
  }
}
Iterator::~Iterator() = default;
void Iterator::RunCluster() const {
  const auto chemical_potential = energy_estimator_.GetChemicalPotential(solvent_element_);
  json clusters_info_array = json::array();
  for (unsigned long long i = 0; i <= final_number_; i += increment_steps_) {
    std::cout << i << " / " << final_number_ << std::endl;
    cfg::Config config;
    if (config_type_ == "config") {
      config = cfg::Config::ReadConfig(std::to_string(i) + ".cfg");
      config.ReassignLatticeVector();
    } else if (config_type_ == "map") {
      config = cfg::Config::ReadMap("lattice.txt",
                                    "element.txt",
                                    "map" + std::to_string(i) + ".txt");
    } else {
      throw std::invalid_argument("Unknown config type: " + config_type_);
    }
    Cluster cluster_finder(config,
                           solvent_element_,
                           element_set_,
                           smallest_cluster_criteria_,
                           solvent_bond_criteria_,
                           energy_estimator_,
                           chemical_potential);

    json clusters_info = json::object();
    clusters_info["index"] = std::to_string(i);
    clusters_info["time"] = filename_time_hashset_.at(i);
    clusters_info["temperature"] = filename_temperature_hashset_.at(i);
    clusters_info["energy"] = filename_energy_hashset_.at(i);
    clusters_info["clusters"] = cluster_finder.GetClustersInfoAndOutput(
        "cluster", std::to_string(i) + "_cluster.cfg");
    clusters_info_array.push_back(clusters_info);

    std::ofstream ofs("cluster_info.json", std::ofstream::out);
    ofs << std::setprecision(16) << clusters_info_array.dump(2) << std::endl;
  }
}
void Iterator::RunShortRangeOrder() const {
  // start
  std::ofstream ofs1("sro1_log.txt", std::ofstream::out);
  std::ofstream ofs2("sro2_log.txt", std::ofstream::out);
  std::ofstream ofs3("sro3_log.txt", std::ofstream::out);

  for (unsigned long long i = 0; i <= final_number_; i += increment_steps_) {
    std::cout << i << " / " << final_number_ << std::endl;
    cfg::Config config;
    if (config_type_ == "config") {
      config = cfg::Config::ReadConfig(std::to_string(i) + ".cfg");
      config.ReassignLatticeVector();
    } else if (config_type_ == "map") {
      config = cfg::Config::ReadMap("lattice.txt",
                                    "element.txt",
                                    "map" + std::to_string(i) + ".txt");
    } else {
      throw std::invalid_argument("Unknown config type: " + config_type_);
    }
    ShortRangeOrder short_range_order(config, element_set_);
    for (size_t j = 1; j <= 3; ++j) {
      std::ofstream *ofs;
      switch (j) {
        case 1:ofs = &ofs1;
          break;
        case 2:ofs = &ofs2;
          break;
        case 3:ofs = &ofs3;
          break;
        default:throw std::invalid_argument("Unknown short range order type: " + std::to_string(j));
      }
      auto sro_map = short_range_order.FindWarrenCowley(j);
      if (i == 0) {
        *ofs << "index\ttime\ttemperature\tenergy\t";
        std::transform(sro_map.cbegin(),
                       sro_map.cend(),
                       std::experimental::make_ostream_joiner(*ofs, "\t"),
                       [](const auto &ii) { return ii.first; });
        *ofs << std::endl;
      }
      *ofs << i << '\t'
           << filename_time_hashset_.at(i) << '\t'
           << filename_temperature_hashset_.at(i) << '\t'
           << filename_energy_hashset_.at(i) << '\t';
      std::transform(sro_map.cbegin(),
                     sro_map.cend(),
                     std::experimental::make_ostream_joiner(*ofs, "\t"),
                     [](const auto &ii) { return ii.second; });
      *ofs << std::endl;
    }
  }
}
void Iterator::RunReformat() const {
  for (unsigned long long i = 0; i <= final_number_; i += increment_steps_) {
    std::cout << i << " / " << final_number_ << std::endl;
    if (config_type_ == "map") {
      auto config = cfg::Config::ReadMap("lattice.txt",
                                         "element.txt",
                                         "map" + std::to_string(i) + ".txt");
      config.WriteConfig(std::to_string(i) + ".cfg", false);
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
} // ansys