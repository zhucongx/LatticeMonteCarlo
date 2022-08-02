#include "Parameter.h"
#include "Utility.h"

namespace api {
Parameter::Parameter() = default;
Parameter::Parameter(int argc, char *argv[]) {
  ParseArgs(argc, argv);
  ReadParam(parameters_filename);
}
Parameter::Parameter(const std::string &param_filename) {
  ReadParam(param_filename);
}
void Parameter::ParseArgs(int argc, char *argv[]) {
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      parameters_filename = std::string(argv[++i]);
  }
}
void Parameter::ReadParam(const std::string &param_filename) {
  if (param_filename.empty()) {
    return;
  }
  std::ifstream ifs(param_filename, std::ifstream::in);
  if (!ifs.is_open()) std::cerr << " error opening " << param_filename << std::endl;
  std::string buffer;
  while (getline(ifs, buffer)) {
    if (buffer.empty()) {
      continue;
    }
    if (buffer[0] == '#') {
      continue;
    }
    std::vector<std::string> segs(split(buffer, " "));
    if (segs[0] == "simulation_method") {
      method = std::string(segs[1]);
    } else if (segs[0] == "config_filename") {
      config_filename_ = std::string(segs[1]);
    } else if (segs[0] == "map_filename") {
      map_filename_ = std::string(segs[1]);
    } else if (segs[0] == "json_coefficients_filename") {
      json_coefficients_filename_ = std::string(segs[1]);
    } else if (segs[0] == "log_dump_steps") {
      log_dump_steps_ = stoull(segs[1]);
    } else if (segs[0] == "config_dump_steps") {
      config_dump_steps_ = stoull(segs[1]);
    } else if (segs[0] == "maximum_steps") {
      maximum_steps_ = stoull(segs[1]);
    } else if (segs[0] == "temperature") {
      temperature_ = stod(segs[1]);
    } else if (segs[0] == "element_set") {
      element_set_.clear();
      std::copy(segs.begin() + 1, segs.end(),
                std::back_inserter(element_set_));
    } else if (segs[0] == "initial_steps") {
      initial_steps_ = stoull(segs[1]);
    } else if (segs[0] == "increment_steps") {
      increment_steps_ = stoull(segs[1]);
    } else if (segs[0] == "smallest_cluster_criteria") {
      smallest_cluster_criteria_ = stoul(segs[1]);
    } else if (segs[0] == "restart_steps") {
      restart_steps_ = stoull(segs[1]);
    } else if (segs[0] == "restart_energy") {
      restart_energy_ = stod(segs[1]);
    } else if (segs[0] == "restart_time") {
      restart_time_ = stod(segs[1]);
    } else if (segs[0] == "factor") {
      factor_ = stoul(segs[1]);
    } else if (segs[0] == "solvent_element") {
      solvent_element_ = std::string(segs[1]);
    } else if (segs[0] == "solute_element_set") {
      solute_element_set_.clear();
      std::copy(segs.begin() + 1, segs.end(),
                std::back_inserter(solute_element_set_));
    } else if (segs[0] == "solute_number_set") {
      solute_number_set_.clear();
      std::transform(segs.begin() + 1,
                     segs.end(),
                     std::back_inserter(solute_number_set_),
                     [](const auto &number) { return stoul(number); });
    } else if (segs[0] == "early_stop_steps") {
      early_stop_steps_ = stoull(segs[1]);
    }
  }
  ifs.close();
}
void Parameter::PrintParameters() const {
  std::cout << "Parameters" << std::endl;
  std::cout << "simulation_method: " << method << std::endl;
  if (method == "FirstKmc" || method == "ChainKmc") {
    std::cout << "config_filename: " << config_filename_ << std::endl;
    std::cout << "map_filename: " << map_filename_ << std::endl;
    std::cout << "log_dump_steps: " << log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << maximum_steps_ << std::endl;
    std::cout << "temperature: " << temperature_ << std::endl;
    std::cout << "element_set: ";
    std::copy(element_set_.begin(), element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "restart_steps: " << restart_steps_ << std::endl;
    std::cout << "restart_energy: " << restart_energy_ << std::endl;
    std::cout << "restart_time: " << restart_time_ << std::endl;
    std::cout << "json_coefficients_filename: " << json_coefficients_filename_ << std::endl;
  } else if (method == "SimulatedAnnealing") {
    std::cout << "factor: " << factor_ << std::endl;
    std::cout << "solvent_element: " << solvent_element_ << std::endl;
    std::cout << "solute_element_set: ";
    std::copy(solute_element_set_.begin(), solute_element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "solute_number_set: ";
    std::transform(solute_number_set_.begin(), solute_number_set_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto number) { return std::to_string(number); });
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << maximum_steps_ << std::endl;
    std::cout << "early_stop_steps: " << early_stop_steps_ << std::endl;
    std::cout << "temperature: " << temperature_ << std::endl;
    std::cout << "json_coefficients_filename: " << json_coefficients_filename_ << std::endl;
  } else if (method == "Cluster") {
    std::cout << "solvent_element: " << solvent_element_ << std::endl;
    std::cout << "element_set: ";
    std::copy(element_set_.begin(), element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "initial_steps: " << initial_steps_ << std::endl;
    std::cout << "increment_steps: " << increment_steps_ << std::endl;
    std::cout << "smallest_cluster_criteria: " << smallest_cluster_criteria_ << std::endl;
    std::cout << "json_coefficients_filename: " << json_coefficients_filename_ << std::endl;
  }
}

} // namespace api
