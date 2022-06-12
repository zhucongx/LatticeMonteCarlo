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
    std::vector<std::string> segs(split(buffer, " "));
    if (segs[0] == "simulation_method") {
      method = segs[1];
    } else if (segs[0] == "config_filename") {
      config_filename_ = segs[1];
    } else if (segs[0] == "json_coefficients_filename") {
      json_coefficients_filename_ = segs[1];
    } else if (segs[0] == "log_dump_steps") {
      log_dump_steps_ = stoull(segs[1]);
    } else if (segs[0] == "config_dump_steps") {
      config_dump_steps_ = stoull(segs[1]);
    } else if (segs[0] == "maximum_number") {
      maximum_number_ = stoull(segs[1]);
    } else if (segs[0] == "temperature") {
      temperature_ = stod(segs[1]);
    } else if (segs[0] == "element_symbols") {
      element_symbols_.clear();
      std::copy(segs.begin() + 1, segs.end(),
                std::back_inserter(element_symbols_));
    } else if (segs[0] == "restart_steps") {
      restart_steps_ = stoull(segs[1]);
    } else if (segs[0] == "restart_energy") {
      restart_energy_ = stod(segs[1]);
    } else if (segs[0] == "restart_time") {
      restart_time_ = stod(segs[1]);
    }
  }
  ifs.close();
}
void Parameter::PrintParameters() const {
  std::cout << "KMC Parameters" << std::endl;
  std::cout << "simulation_method: " << method << std::endl;
  std::cout << "config_filename: " << config_filename_ << std::endl;
  std::cout << "json_coefficients_filename: " << json_coefficients_filename_ << std::endl;
  std::cout << "log_dump_steps: " << log_dump_steps_ << std::endl;
  std::cout << "config_dump_steps: " << config_dump_steps_ << std::endl;
  std::cout << "maximum_number: " << maximum_number_ << std::endl;
  std::cout << "temperature: " << temperature_ << std::endl;

  std::cout << "element_symbols: ";
  std::copy(element_symbols_.begin(), element_symbols_.end(),
            std::ostream_iterator<std::string>(std::cout, " "));
  std::cout << std::endl;
  std::cout << "restart_steps: " << restart_steps_ << std::endl;
  std::cout << "restart_energy: " << restart_energy_ << std::endl;
  std::cout << "restart_time: " << restart_time_ << std::endl;
}

} // namespace api
