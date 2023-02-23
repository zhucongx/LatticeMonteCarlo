#include "Iterator.h"

#include <utility>
#include <algorithm>
#include <iostream>
#include <omp.h>
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
#pragma omp parallel  default(none) shared(std::cout)
  {
#pragma omp master
    {
      std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    }
  }
}
Iterator::~Iterator() = default;
void Iterator::RunAnsys() const {
  const auto chemical_potential = energy_estimator_.GetChemicalPotential(solvent_element_);
  json ansys_info_array = json::array();
#pragma omp parallel for default(none) shared(ansys_info_array, chemical_potential, std::cout)
  for (unsigned long long i = initial_steps_; i <= final_number_; i += increment_steps_) {
#pragma omp critical
    {
      std::cout << i << " / " << final_number_ << std::endl;
    }
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
    // basic information
    json ansys_info = json::object();
    ansys_info["index"] = std::to_string(i);
    ansys_info["time"] = filename_time_hashset_.at(i);
    ansys_info["temperature"] = filename_temperature_hashset_.at(i);
    ansys_info["energy"] = filename_energy_hashset_.at(i);
    // cluster information
    ansys_info["clusters"] =
        Cluster(config,
                solvent_element_,
                element_set_,
                smallest_cluster_criteria_,
                solvent_bond_criteria_,
                energy_estimator_,
                chemical_potential).GetClustersInfoAndOutput(
            "cluster", std::to_string(i) + "_cluster.cfg");
    // sro information
    ShortRangeOrder short_range_order(config, element_set_);
    ansys_info["short_range_order"]["first"] = short_range_order.FindWarrenCowley(1);
    ansys_info["short_range_order"]["second"] = short_range_order.FindWarrenCowley(2);
    ansys_info["short_range_order"]["third"] = short_range_order.FindWarrenCowley(3);
#pragma omp critical
    {
      ansys_info_array.push_back(ansys_info);
    }
  }
  std::ofstream ofs("ansys_info.json", std::ofstream::out);
  ofs.precision(16);
  ofs << ansys_info_array.dump(2) << std::endl;
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