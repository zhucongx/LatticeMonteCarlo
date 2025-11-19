#include "Traverse.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
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
                   double escape_temperature,
                   const std::string &predictor_filename,
                   std::string log_type,
                   std::string config_type)
    : initial_steps_(initial_steps),
      increment_steps_(increment_steps),
      final_steps_(increment_steps),
      solvent_element_(solvent_element),
      element_set_(std::move(element_set)),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria),
      escape_temperature_(escape_temperature),
      log_map_{},
      log_type_(std::move(log_type)),
      config_type_(std::move(config_type)),
      energy_predictor_(predictor_filename, element_set_),
      vacancy_migration_predictor_(predictor_filename, GetConfig(config_type_, initial_steps_), element_set_),
      energy_change_predictor_pair_site_(predictor_filename, GetConfig(config_type_, initial_steps_), element_set_) {
  std::string log_file_name;
  if (log_type_ == "kinetic_mc") {
    log_file_name = "kmc_log.txt";
  } else if (log_type_ == "canonical_mc") {
    log_file_name = "cmc_log.txt";
  } else if (log_type_ == "simulated_annealing") {
    log_file_name = "sa_log.txt";
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
    final_steps_ = step_number;
    size_t col_index = 1;
    while (line_stream >> buffer) {
      const auto &key = headers[col_index];
      try {
        const auto double_value = boost::lexical_cast<double>(buffer);
        if (!std::holds_alternative<std::unordered_map<unsigned long long, double>>(log_map_[key])) {
          log_map_[key] = std::unordered_map<unsigned long long, double>();
        }
        std::get<std::unordered_map<unsigned long long, double>>(log_map_[key])[step_number] = double_value;
      } catch (const boost::bad_lexical_cast &) {
        if (!std::holds_alternative<std::unordered_map<unsigned long long, std::string>>(log_map_[key])) {
          log_map_[key] = std::unordered_map<unsigned long long, std::string>();
        }
        std::get<std::unordered_map<unsigned long long, std::string>>(log_map_[key])[step_number] = buffer;
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
  std::ofstream  frame_ofs("ansys_frame_log.txt", std::ofstream::out);
  std::ofstream  cluster_ofs("ansys_cluster_log.txt", std::ofstream::out);
  const auto chemical_potential = energy_predictor_.GetChemicalPotential(solvent_element_);
  // nlohmann::json ansys_info_array = nlohmann::json::array();
  frame_ofs << GetHeaderFrameString() << std::flush;
  cluster_ofs << GetHeaderClusterString() << std::flush;
#pragma omp parallel for default(none) schedule(static, 1) shared(chemical_potential, std::cout, std::cerr, frame_ofs, cluster_ofs)
  for (unsigned long long i = initial_steps_; i <= final_steps_; i += increment_steps_) {
#pragma omp critical
    {
      std::cout << i << " / " << final_steps_ << " " << std::fixed << std::setprecision(2)
                << static_cast<double>(i - initial_steps_ + 1) /
              static_cast<double>(final_steps_ - initial_steps_ + 1) * 100
                << "%" << std::endl;
    }
    cfg::Config config = GetConfig(config_type_, i);

    // basic information
    nlohmann::json frame_info = nlohmann::json::object();
    std::map<std::string, cfg::Config::VectorVariant> auxiliary_lists{};
    std::map<std::string, cfg::Config::ValueVariant> global_list;

    frame_info["steps"] = i;
    global_list["steps"] = i;
    const auto time = log_map_.find("time") == log_map_.end()
        ? nan("")
        : std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("time")).at(i);
    frame_info["time"] = time;
    global_list["time"] = time;
    const auto current_temperature = std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("temperature")).at(i);
    frame_info["temperature"] = current_temperature;
    global_list["temperature"] = current_temperature;
    const auto energy = std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("energy")).at(i);
    frame_info["energy"] = energy;
    global_list["energy"] = energy;

    // basic cluster information
    SoluteCluster(config,
                  solvent_element_,
                  element_set_,
                  smallest_cluster_criteria_,
                  solvent_bond_criteria_,
                  energy_predictor_,
                  chemical_potential)
        .GetClustersInfo(frame_info, auxiliary_lists, global_list);

    // binding energy
    ExitTime(config,
             solvent_element_,
             element_set_,
             escape_temperature_,
             current_temperature,
             energy_predictor_,
             vacancy_migration_predictor_,
             energy_change_predictor_pair_site_,
             chemical_potential)
        .GetExitTimeInfo(frame_info, auxiliary_lists, global_list);

    // sro information
    // ShortRangeOrder short_range_order(config, element_set_);
    // const auto sro1 = short_range_order.FindWarrenCowley(1);
    // frame_info["warren_cowley"]["first"] = sro1;
    // for (const auto &[pair, value]: sro1) {
    //   global_list["sro1_" + pair] = value;
    // }
    //
    // const auto sro2 = short_range_order.FindWarrenCowley(2);
    // frame_info["warren_cowley"]["second"] = sro2;
    // for (const auto &[pair, value]: sro2) {
    //   global_list["sro2_" + pair] = value;
    // }
    //
    // const auto sro3 = short_range_order.FindWarrenCowley(3);
    // frame_info["warren_cowley"]["third"] = sro3;
    // for (const auto &[pair, value]: sro3) {
    //   global_list["sro3_" + pair] = value;
    // }

    // vacancy information
    // auto vacancy_lattice_id = config.GetVacancyLatticeId();
    // static auto convert = [](const std::set<Element> &element_set_base,
    //                          const std::map<Element, size_t> &element_map) -> std::map<std::string, size_t> {
    //   std::map<std::string, size_t> string_map;
    //   for (const auto &element: element_set_base) {
    //     string_map[element.GetString()] = 0;
    //   }
    //   for (const auto &[element, value]: element_map) {
    //     string_map[element.GetString()] = value;
    //   }
    //   return string_map;
    // };
    // const auto vac1 = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 1));
    // frame_info["vacancy_local"]["first"] = vac1;
    // for (const auto &[element, value]: vac1) {
    //   global_list["vac1_" + element] = value;
    // }
    // const auto vac2 = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 2));
    // frame_info["vacancy_local"]["second"] = vac2;
    // for (const auto &[element, value]: vac2) {
    //   global_list["vac2_" + element] = value;
    // }
    // const auto vac3 = convert(element_set_, config.GetLocalInfoOfLatticeId(vacancy_lattice_id, 3));
    // frame_info["vacancy_local"]["third"] = vac3;
    // for (const auto &[element, value]: vac3) {
    //   global_list["vac3_" + element] = value;
    // }


    boost::filesystem::create_directories("ansys");
    auxiliary_lists["map_shift"] = config.GetMapShiftList();
    config.WriteExtendedXyz("ansys/" + std::to_string(i) + ".xyz.gz", auxiliary_lists, global_list);
    const auto cluster_string = GetClusterString(frame_info);
    const auto frame_string = GetFrameString(frame_info);
#pragma omp critical
    {
      // ansys_info_array.push_back(frame_info);
      cluster_ofs << cluster_string << std::flush;
      frame_ofs << frame_string << std::flush;
    }
  }
  // std::cout << "Finished. Sorting..." << std::endl;
  // std::sort(ansys_info_array.begin(), ansys_info_array.end(), [](const nlohmann::json &lhs, const nlohmann::json &rhs) {
  //   return lhs["steps"] < rhs["steps"];
  // });
  // std::ofstream ofs("ansys_info.json.gz", std::ios_base::out | std::ios_base::binary);
  // boost::iostreams::filtering_ostream fos;
  // fos.push(boost::iostreams::gzip_compressor());
  // fos.push(ofs);
  // fos.precision(16);
  // fos << ansys_info_array.dump(2) << std::endl;
  std::cout << "Done..." << std::endl;
}

std::string Traverse::GetHeaderClusterString() const {
  std::string header_frame = "steps\ttime\ttemperature\tenergy\tcluster_id\tcluster_energy\tcluster_size\t";
  for (const auto &element: element_set_) {
    header_frame += "cluster_";
    header_frame += element.GetString();
    header_frame += "\t";
  }
  header_frame += "cluster_X\teffective_radius\tmass_gyration_radius\tasphericity\tacylindricity\tanisotropy\t"
                  "markov_escape_time\t"
                  // "markov_current_time\t"
                  "effective_binding_energy\t"
                  "vacancy_profile_energy\tvacancy_binding_energy\t";
  for (const auto &element: element_set_) {
    if (element == Element("X")) {
      continue;
    }
    header_frame += "vacancy_binding_energy_";
    header_frame += element.GetString();
    header_frame += "\t";
  }
  header_frame += "barrier_mean_in\tbarrier_mean_to\tbarrier_mean_on\tbarrier_mean_off\t"
                  "barrier_std_in\tbarrier_std_to\tbarrier_std_on\tbarrier_std_off\t"
                  "geometry_center\n";
  return header_frame;
}

std::string Traverse::GetClusterString(const nlohmann::json &frame) const {
  std::stringstream cluster_stream;
  cluster_stream.precision(16);
  for (const auto &cluster: frame["clusters"]) {
    cluster_stream << frame["steps"] << "\t" << frame["time"] << "\t" << frame["temperature"] << "\t" << frame["energy"]
                   << "\t" << cluster["cluster_id"] << "\t" << cluster["cluster_energy"] << "\t"
                   << cluster["cluster_size"] << "\t";
    for (const auto &element: element_set_) {
      cluster_stream << cluster["elements_number"][element.GetString()] << "\t";
    }
    cluster_stream << cluster["elements_number"]["X"] << "\t" << cluster["effective_radius"] << "\t"
                   << cluster["mass_gyration_radius"] << "\t" << cluster["asphericity"] << "\t"
                   << cluster["acylindricity"] << "\t" << cluster["anisotropy"] << "\t"
                   << cluster["markov_escape_time"] << "\t"
                   // << cluster["markov_current_time"] << "\t"
                   << cluster["effective_binding_energy"] << "\t"
                   << cluster["vacancy_profile_energy"] << "\t" << cluster["vacancy_binding_energy"] << "\t";
    for (const auto &element: element_set_) {
      cluster_stream << cluster["vacancy_binding_energy_" + element.GetString()] << "\t";
    }
    const auto barriers = cluster["barriers"].get<std::vector<double>>();
    cluster_stream << barriers[0] << "\t" << barriers[1] << "\t" << barriers[2] << "\t" << barriers[3] << "\t";
    cluster_stream << barriers[4] << "\t" << barriers[5] << "\t" << barriers[6] << "\t" << barriers[7] << "\t";
    cluster_stream << "[" << GetVectorToString(cluster["geometry_center"], ",") << "]\n";
  }
  return cluster_stream.str();
}

std::string Traverse::GetHeaderFrameString() const {
  std::string header_frame = "steps\ttime\ttemperature\tenergy\tnum_cluster\tnum_atom\t";
  for (const auto &element: element_set_) {
    header_frame += "num_" + element.GetString() + "\t";
  }
  header_frame += "profile_vacancy_concentration\tbinding_vacancy_concentration\tprofile_simp_vacancy_concentration\t"
                  "binding_simp_vacancy_concentration\teffective_binding_simp_vacancy_concentration\t";
  header_frame += "cluster_size_list\t";
  // static const std::vector<std::string> order_list{"first", "second", "third"};
  // for (const auto &order: order_list) {
  //   for (auto element1: element_set_) {
  //     for (auto element2: element_set_) {
  //       header_frame += "warren_cowley_" + order + "_" + element1.GetString() + "-" + element2.GetString() + "\t";
  //     }
  //   }
  // }
  // for (const auto &order: order_list) {
  //   for (const auto &element: element_set_) {
  //     header_frame += "vacancy_local_" + order + "_" + element.GetString() + "\t";
  //   }
  // }
  // for (const auto &element: element_set_) {
  //   header_frame += "vacancy_local_binding_energy_" + element.GetString() + "\t";
  // }
  if (!header_frame.empty() && header_frame.back() == '\t') {
    header_frame.back() = '\n';
  }
  return header_frame;
}

std::string Traverse::GetFrameString(const nlohmann::json &frame) const {
  std::stringstream frame_stream;
  frame_stream.precision(16);

  frame_stream << frame["steps"] << "\t" << frame["time"] << "\t" << frame["temperature"] << "\t" << frame["energy"]
               << "\t" << frame["num_cluster"] << "\t" << frame["num_atom"] << "\t";
  for (const auto &element: element_set_) {
    frame_stream << frame["num_" + element.GetString()] << "\t";
  }
  frame_stream << frame["profile_vacancy_concentration"] << "\t"
               << frame["binding_vacancy_concentration"] << "\t"
               << frame["profile_simp_vacancy_concentration"] << "\t"
               << frame["binding_simp_vacancy_concentration"] << "\t"
               << frame["effective_binding_simp_vacancy_concentration"] << "\t";
  frame_stream << "[" << boost::algorithm::join(frame["cluster_size_list"].get<std::vector<std::string>>(), ",")
               << "]\t";

  // for (const auto &order: frame["warren_cowley"]) {
  //   for (const auto &pair: order.items()) {
  //     frame_stream << pair.value() << "\t";
  //   }
  // }
  // for (const auto &order: frame["vacancy_local"]) {
  //   for (const auto &pair: order.items()) {
  //     frame_stream << pair.value() << "\t";
  //   }
  // }
  // for (const auto &element: element_set_) {
  //   frame_stream << frame["vacancy_local_binding_energy"][element.GetString()] << "\t";
  // }

  auto frame_string = frame_stream.str();
  if (!frame_string.empty() && frame_string.back() == '\t') {
    frame_string.back() = '\n';
  }
  return frame_string;
}

void Traverse::RunReformat() const {
#pragma omp parallel for default(none) schedule(static, 1) shared(std::cout, std::cerr)
  for (unsigned long long i = initial_steps_; i <= final_steps_; i += increment_steps_) {
#pragma omp critical
    {
      std::cout << i << " / " << final_steps_ << " " << std::fixed << std::setprecision(2)
                << static_cast<double>(i - initial_steps_ + 1) /
              static_cast<double>(final_steps_ - initial_steps_ + 1) * 100
                << "%" << std::endl;
    }
    // try {
      if (config_type_ == "map") {
        auto config = cfg::Config::ReadMap("lattice.txt", "element.txt", "map" + std::to_string(i) + ".txt");
        config.WriteConfig(std::to_string(i) + ".cfg.gz");
      } else if (config_type_ == "config") {
        auto config = cfg::Config::ReadConfig(std::to_string(i) + ".cfg.gz");
        config.ReassignLatticeVector();
        if (i == 0) {
          config.WriteLattice("lattice.txt");
          config.WriteElement("element.txt");
        }
        config.WriteMap("map" + std::to_string(i) + ".txt");
      } else {
        throw std::invalid_argument("Unknown config type: " + config_type_);
      }
    // } catch (const std::exception &e) {
    //   std::cerr << "Error when working with " + std::to_string(i) << std::endl;
    // }
  }
  std::cout << "Done..." << std::endl;
}
}    // namespace ansys