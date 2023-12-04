/**************************************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 6/4/22 4:53 AM                                                                          *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/27/23 3:30 PM                                                           *
 **************************************************************************************************/

#include "Traverse.h"

#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace ansys {

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
  if (config_type_ != "config" && config_type_ != "map") {
    throw std::invalid_argument("Unknown config type: " + config_type_);
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
#pragma omp parallel default(none) shared(std::cout)
  {
#pragma omp master
    {
      std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    }
  }
}

Traverse::~Traverse() = default;

void Traverse::RunAnsys() const {
  const auto chemical_potential = energy_estimator_.GetChemicalPotential(solvent_element_);
  json ansys_info_array = json::array();
#pragma omp parallel default(none) shared(ansys_info_array, chemical_potential, std::cout)
  {
#pragma omp for schedule(static, 1)
    for (unsigned long long i = initial_steps_; i <= final_number_; i += increment_steps_) {
#pragma omp critical
      {
        std::cout << i << " / " << final_number_ << std::endl;
      }
      Config config;
      if (config_type_ == "config") {
        config = Config::ReadCfg(std::to_string(i) + ".cfg");
      } else if (config_type_ == "map") {
        config = Config::ReadMap("lattice.txt",
                                 "element.txt",
                                 "map" + std::to_string(i) + ".txt");
      }
      // basic information
      json ansys_info = json::object();
      auto time = filename_time_hashset_.at(i);
      auto temperature = filename_temperature_hashset_.at(i);
      auto energy = filename_energy_hashset_.at(i);
      ansys_info["index"] = i;
      ansys_info["time"] = time;
      ansys_info["temperature"] = temperature;
      ansys_info["energy"] = energy;

      std::map<std::string, Config::ValueVariant> global_info_map{
          {"index", i}, {"time", time}, {"temperature", temperature}, {"energy", energy}};
      // cluster information
      ansys_info["clusters"] =
          SoluteCluster(config,
                        solvent_element_,
                        element_set_,
                        smallest_cluster_criteria_,
                        solvent_bond_criteria_,
                        energy_estimator_,
                        chemical_potential).GetClustersInfoAndOutput(
              "cluster", std::to_string(i) + "_cluster.cfg", global_info_map);
      // // sro information
      // ShortRangeOrder short_range_order(config, element_set_);
      // ansys_info["warren_cowley"]["first"] = short_range_order.FindWarrenCowley(1);
      // ansys_info["warren_cowley"]["second"] = short_range_order.FindWarrenCowley(2);
      // ansys_info["warren_cowley"]["third"] = short_range_order.FindWarrenCowley(3);
#pragma omp critical
      {
        ansys_info_array.push_back(ansys_info);
      }
      if (omp_get_thread_num() == 0 && omp_get_num_threads() == 1) {

        std::ofstream ofs("ansys_info.json.gz", std::ios_base::out | std::ios_base::binary);
        boost::iostreams::filtering_ostream fos;
        fos.push(boost::iostreams::gzip_compressor());
        fos.push(ofs);
        fos.precision(16);
        fos << ansys_info_array.dump(2) << std::endl;
      }
    }
  }
  std::cout << "Finished. Sorting..." << std::endl;
  std::sort(ansys_info_array.begin(), ansys_info_array.end(),
            [](const json &lhs, const json &rhs) {
              return lhs["index"] < rhs["index"];
            });
  std::ofstream ofs("ansys_info.json", std::ofstream::out);
  ofs.precision(16);
  ofs << ansys_info_array.dump(2) << std::endl;
  std::cout << "Done..." << std::endl;
}
void Traverse::RunReformat() const {
#pragma omp parallel default(none) shared(std::cout)
  {
#pragma omp for schedule(static, 1)
    for (unsigned long long i = initial_steps_; i <= final_number_; i += increment_steps_) {
      std::cout << std::to_string(i) + " / " + std::to_string(final_number_) << std::endl;
      if (config_type_ == "map") {
        // Todo: catch exception where map file does not exist and break the loop
        auto config = Config::ReadMap("lattice.txt", "element.txt", "map" + std::to_string(i) + ".txt");
        config.WriteConfig(std::to_string(i) + ".cfg.gz");
      } else {
        throw std::invalid_argument("Unknown config type: " + config_type_);
      }
    }
  }
  std::cout << "Done..." << std::endl;
}
} // ansys