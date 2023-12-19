#include "Traverse.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <omp.h>

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
      log_json_(nlohmann::json::object()),
      log_type_(std::move(log_type)),
      config_type_(std::move(config_type)),
      energy_predictor_(predictor_filename, element_set_),
      vacancy_migration_predictor_(predictor_filename, GetConfig(config_type_, 0), element_set_),
      energy_change_predictor_pair_site_(predictor_filename, GetConfig(config_type_, 0), element_set_),
      frame_ofs_("ansys_frame_log.txt", std::ofstream::out),
      cluster_ofs_("ansys_cluster_log.txt", std::ofstream::out) {
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
  while (ifs.peek() == '#') {
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  std::string buffer;
  // read header
  std::getline(ifs, buffer);
  std::vector<std::string> headers;
  boost::algorithm::split(headers, buffer, boost::is_any_of("\t"));
  for (const auto &header: headers) {
    if (header=="steps") {
      continue;
    }
    log_json_[header] = nlohmann::json::object();
  }
  // read data
  while (std::getline(ifs, buffer)) {
    if (buffer.empty()) {
      continue;
    }
    if (buffer[0] == '#') {
      continue;
    }
    std::istringstream line_stream(buffer);
    unsigned long long step_number;
    line_stream >> step_number;
    if (step_number < initial_steps_ || (step_number - initial_steps_) % increment_steps != 0) {
      continue;
    }
    final_number_ = step_number;
    size_t col_index = 1;
    while (line_stream >> buffer) {
      try {
        const auto double_value = boost::lexical_cast<double>(buffer);
        log_json_[headers[col_index]][step_number] = double_value;
      } catch (const boost::bad_lexical_cast &) {
        log_json_[headers[col_index]][step_number] = buffer;
      }
      col_index++;
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
  nlohmann::json ansys_info_array = nlohmann::json::array();
  frame_ofs_ << GetHeaderFrameString() << std::flush;
  cluster_ofs_ << GetHeaderClusterString() << std::flush;
#pragma omp parallel for default(none) schedule(static, 1) shared(ansys_info_array, chemical_potential, std::cout)
  for (unsigned long long i = initial_steps_; i <= final_number_; i += increment_steps_) {
#pragma omp critical
    {
      std::cout << i << " / " << final_number_ << " " << std::fixed << std::setprecision(2)
                << static_cast<double>(i) / static_cast<double>(final_number_) * 100 << "%" << std::endl;
    }
    cfg::Config config = GetConfig(config_type_, i);

    // basic information
    nlohmann::json ansys_info = nlohmann::json::object();
    std::map<std::string, cfg::Config::ValueVariant> global_list;
    ansys_info["steps"] = i;
    global_list["steps"] = i;

    ansys_info["time"] = log_json_["time"].at(i);
    global_list["time"] = log_json_["time"].at(i).get<double>();

    ansys_info["temperature"] = log_json_["temperature"].at(i);
    global_list["temperature"] = log_json_["temperature"].at(i).get<double>();

    ansys_info["energy"] = log_json_["energy"].at(i);
    global_list["energy"] = log_json_["energy"].at(i).get<double>();

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
    ansys_info["vacancy_local"]["first"] = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 1));
    ansys_info["vacancy_local"]["second"] =
        convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 2));
    ansys_info["vacancy_local"]["third"] = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 3));

    // binding energy
    const auto exit_time = ExitTime(config,
                                    solvent_element_,
                                    log_json_["temperature"].at(i),
                                    vacancy_migration_predictor_,
                                    energy_change_predictor_pair_site_,
                                    chemical_potential);
    const auto binding_energy = exit_time.GetBindingEnergy();
    auxiliary_lists["binding_energy"] = binding_energy;

    ansys_info["vacancy_local_binding_energy"] = binding_energy[config.GetVacancyAtomId()];
    global_list["vacancy_local_binding_energy"] = binding_energy[config.GetVacancyAtomId()];

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
    const auto cluster_string = GetClusterString(ansys_info);
    const auto frame_string = GetFrameString(ansys_info);
#pragma omp critical
    {
      ansys_info_array.push_back(ansys_info);
      cluster_ofs_ << cluster_string << std::flush;
      frame_ofs_ << frame_string << std::flush;
    }
  }
  std::cout << "Finished. Sorting..." << std::endl;
  std::sort(ansys_info_array.begin(), ansys_info_array.end(), [](const nlohmann::json &lhs, const nlohmann::json &rhs) {
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

std::string Traverse::GetHeaderClusterString() const {
  std::string header_frame = "steps\ttime\ttemperature\tenergy\tcluster_energy\tcluster_size\t";
  for (const auto &element: element_set_) {
    header_frame += "cluster_";
    header_frame += element.GetString();
    header_frame += "\t";
  }
  header_frame += "cluster_X\tmass_gyration_radius\tasphericity\tacylindricity\tanisotropy\t"
                  "vacancy_binding_energy\tvacancy_local_binding_energy\n";
  return header_frame;
}

std::string Traverse::GetClusterString(const nlohmann::json &frame) const {
  std::stringstream cluster_stream;
  cluster_stream.precision(16);
  for (const auto &cluster: frame["clusters"]) {
    cluster_stream << frame["steps"] << "\t" << frame["time"] << "\t" << frame["temperature"] << "\t" << frame["energy"]
                   << "\t" << cluster["cluster_energy"] << "\t" << cluster["cluster_size"] << "\t";
    for (const auto &element: element_set_) {
      cluster_stream << cluster["elements_number"][element.GetString()] << "\t";
    }
    cluster_stream << cluster["elements_number"]["X"] << "\t" << cluster["mass_gyration_radius"] << "\t"
                   << cluster["shape"]["asphericity"] << "\t" << cluster["shape"]["acylindricity"] << "\t"
                   << cluster["shape"]["anisotropy"] << "\t" << cluster["vacancy_binding_energy"] << "\t"
                   << frame["vacancy_local_binding_energy"] << "\n";
  }
  return cluster_stream.str();
}

std::string Traverse::GetHeaderFrameString() const {
  std::string header_frame = "steps\ttime\ttemperature\tenergy\tnum_cluster\tnum_atom\t";
  for (const auto &element: element_set_) {
    header_frame += "num_" + element.GetString() + "\t";
  }
  header_frame += "cluster_size_list\t";
  static const std::vector<std::string> order_list{"first", "second", "third"};
  for (const auto &order: order_list) {
    for (auto it1 = element_set_.cbegin(); it1 != element_set_.cend(); ++it1) {
      for (auto it2 = element_set_.cbegin(); it2 != element_set_.cend(); ++it2) {
        header_frame += "warren_cowley_" + order + "_" + it1->GetString() + "-" + it2->GetString() + "\t";
      }
    }
  }
  for (const auto &order: order_list) {
    for (const auto &element: element_set_) {
      header_frame += "vacancy_local_" + order + "_" + element.GetString() + "\t";
    }
  }
  header_frame += "vacancy_local_binding_energy\n";
  return header_frame;
}

std::string Traverse::GetFrameString(const nlohmann::json &frame) const {
  std::stringstream frame_stream;
  frame_stream.precision(16);

  frame_stream << frame["steps"] << "\t" << frame["time"] << "\t" << frame["temperature"] << "\t" << frame["energy"]
               << "\t";

  constexpr size_t kCriticalSize = 10;
  size_t num_cluster = 0;
  size_t num_atom = 0;
  std::vector<std::string> cluster_size_list;
  std::map<Element, size_t> element_number_map{};
  for (const auto &element: element_set_) {
    element_number_map[element] = 0;
  }

  for (const auto &cluster: frame["clusters"]) {
    if (cluster["cluster_size"] >= kCriticalSize) {
      num_cluster++;
      num_atom += cluster["cluster_size"].get<size_t>();
      cluster_size_list.push_back(std::to_string(cluster["cluster_size"].get<size_t>()));
      for (const auto &element: element_set_) {
        element_number_map[element] += cluster["elements_number"][element.GetString()].get<size_t>();
      }
    }
  }
  frame_stream << num_cluster << "\t" << num_atom << "\t";
  for (const auto &element: element_set_) {
    frame_stream << element_number_map[element] << "\t";
  }
  frame_stream << "\"[" << boost::algorithm::join(cluster_size_list, ",") << "]\"\t";

  for (const auto &order: frame["warren_cowley"]) {
    for (const auto &pair: order.items()) {
      frame_stream << pair.value() << "\t";
    }
  }
  for (const auto &order: frame["vacancy_local"]) {
    for (const auto &pair: order.items()) {
      frame_stream << pair.value() << "\t";
    }
  }
  frame_stream << frame["vacancy_local_binding_energy"] << "\n";

  return frame_stream.str();
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