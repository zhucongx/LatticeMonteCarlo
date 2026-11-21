#include "Parameter.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace api {
namespace {
unsigned long long int GetTimeBasedSeed() {
  return static_cast<unsigned long long int>(std::chrono::system_clock::now().time_since_epoch().count());
}

// `location` is the "file:line: " prefix of the message.
[[noreturn]] void ThrowInvalidValue(const std::string &location,
                                    const std::string &key,
                                    const std::string &value,
                                    const std::string &expected) {
  throw std::invalid_argument(location + "invalid value '" + value + "' for parameter " + key + ", expected " +
                              expected);
}

// std::stoull alone would wrap "-5" around to a huge number and read "1e9" as 1.
unsigned long long int ParseUnsigned(const std::string &location, const std::string &key, const std::string &value) {
  const bool only_digits = std::all_of(value.begin(), value.end(), [](const unsigned char character) {
    return std::isdigit(character) != 0;
  });
  if (value.empty() || !only_digits) {
    ThrowInvalidValue(location, key, value, "a non-negative integer written with digits only");
  }
  try {
    return std::stoull(value);
  } catch (const std::out_of_range &) {
    ThrowInvalidValue(location, key, value, "an integer that fits in 64 bits");
  }
}

// std::stod alone would read "450K" as 450.
double ParseDouble(const std::string &location, const std::string &key, const std::string &value) {
  size_t num_parsed_characters = 0;
  double number = 0;
  try {
    number = std::stod(value, &num_parsed_characters);
  } catch (const std::exception &) {
    num_parsed_characters = 0;
  }
  if (num_parsed_characters == 0 || num_parsed_characters != value.size() || !std::isfinite(number)) {
    ThrowInvalidValue(location, key, value, "a number");
  }
  return number;
}

bool ParseBool(const std::string &location, const std::string &key, const std::string &value) {
  std::string lower_case(value);
  std::transform(lower_case.begin(), lower_case.end(), lower_case.begin(), [](const unsigned char character) {
    return static_cast<char>(std::tolower(character));
  });
  if (lower_case == "true" || lower_case == "1") {
    return true;
  }
  if (lower_case == "false" || lower_case == "0") {
    return false;
  }
  ThrowInvalidValue(location, key, value, "true or false");
}
}    // namespace

Parameter::Parameter(int argc, char *argv[])
    : seed_(GetTimeBasedSeed()) {
  ParseArgs(argc, argv);
  ReadParam(parameters_filename);
}

Parameter::Parameter(const std::string &param_filename)
    : seed_(GetTimeBasedSeed()) {
  ReadParam(param_filename);
}

void Parameter::ParseArgs(int argc, char *argv[]) {
  for (int i = 0; i < argc; i++) {
    if (!std::strcmp(argv[i], "--p") || !std::strcmp(argv[i], "-p")) parameters_filename = std::string(argv[++i]);
  }
}

void Parameter::ReadParam(const std::string &param_filename) {
  if (param_filename.empty()) {
    return;
  }
  std::ifstream ifs(param_filename, std::ifstream::in);
  if (!ifs.is_open()) {
    throw std::runtime_error("Cannot open " + param_filename);
  }
  std::string buffer;
  size_t line_number = 0;
  while (std::getline(ifs, buffer)) {
    ++line_number;
    // Split on any whitespace (blanks, tabs and the carriage return of a Windows line ending) and stop at a
    // comment. A blank line or a comment line leaves no token.
    std::istringstream line_stream(buffer);
    std::vector<std::string> segs;
    for (std::string token; line_stream >> token;) {
      if (token[0] == '#') {
        break;
      }
      segs.push_back(token);
    }
    if (segs.empty()) {
      continue;
    }
    const std::string location = param_filename + ":" + std::to_string(line_number) + ": ";
    if (segs.size() < 2) {
      throw std::invalid_argument(location + "parameter " + segs[0] + " has no value");
    }
    if (segs[0] == "simulation_method") {
      method = std::string(segs[1]);
    } else if (segs[0] == "config_filename") {
      config_filename_ = std::string(segs[1]);
    } else if (segs[0] == "map_filename") {
      map_filename_ = std::string(segs[1]);
    } else if (segs[0] == "json_coefficients_filename") {
      json_coefficients_filename_ = std::string(segs[1]);
    } else if (segs[0] == "time_temperature_filename") {
      time_temperature_filename_ = std::string(segs[1]);
    } else if (segs[0] == "log_type") {
      log_type_ = std::string(segs[1]);
    } else if (segs[0] == "config_type") {
      config_type_ = std::string(segs[1]);
    } else if (segs[0] == "log_dump_steps") {
      log_dump_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "config_dump_steps") {
      config_dump_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "maximum_steps") {
      maximum_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "thermodynamic_averaging_steps") {
      thermodynamic_averaging_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "temperature") {
      temperature_ = ParseDouble(location, segs[0], segs[1]);
    } else if (segs[0] == "initial_temperature") {
      initial_temperature_ = ParseDouble(location, segs[0], segs[1]);
    } else if (segs[0] == "decrement_temperature") {
      decrement_temperature_ = ParseDouble(location, segs[0], segs[1]);
    } else if (segs[0] == "element_set") {
      element_set_.clear();
      std::copy(segs.begin() + 1, segs.end(), std::back_inserter(element_set_));
    } else if (segs[0] == "initial_steps") {
      initial_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "increment_steps") {
      increment_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "smallest_cluster_criteria") {
      smallest_cluster_criteria_ = static_cast<size_t>(ParseUnsigned(location, segs[0], segs[1]));
    } else if (segs[0] == "solvent_bond_criteria") {
      solvent_bond_criteria_ = static_cast<size_t>(ParseUnsigned(location, segs[0], segs[1]));
    } else if (segs[0] == "escape_temperature") {
      escape_temperature_ = ParseDouble(location, segs[0], segs[1]);
    } else if (segs[0] == "restart_steps") {
      restart_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "rate_corrector") {
      rate_corrector_ = ParseBool(location, segs[0], segs[1]);
    } else if (segs[0] == "early_stop") {
      early_stop_ = ParseBool(location, segs[0], segs[1]);
    } else if (segs[0] == "solute_disp") {
      solute_disp_ = ParseBool(location, segs[0], segs[1]);
    } else if (segs[0] == "restart_energy") {
      restart_energy_ = ParseDouble(location, segs[0], segs[1]);
    } else if (segs[0] == "restart_time") {
      restart_time_ = ParseDouble(location, segs[0], segs[1]);
    } else if (segs[0] == "factor") {
      factor_ = static_cast<size_t>(ParseUnsigned(location, segs[0], segs[1]));
    } else if (segs[0] == "solvent_element") {
      solvent_element_ = std::string(segs[1]);
    } else if (segs[0] == "solute_element_set") {
      solute_element_set_.clear();
      std::copy(segs.begin() + 1, segs.end(), std::back_inserter(solute_element_set_));
    } else if (segs[0] == "solute_number_set") {
      solute_number_set_.clear();
      std::transform(segs.begin() + 1, segs.end(), std::back_inserter(solute_number_set_), [&](const auto &number) {
        return static_cast<size_t>(ParseUnsigned(location, segs[0], number));
      });
    } else if (segs[0] == "early_stop_steps") {
      early_stop_steps_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "seed") {
      seed_ = ParseUnsigned(location, segs[0], segs[1]);
    } else if (segs[0] == "vacancy_trajectory") {
      // Read by earlier versions. Accepted so that their parameter files still load.
      warnings_.push_back(location + "parameter vacancy_trajectory is no longer used and is ignored");
    } else {
      // A misspelled key must not fall back to the default value of the parameter it was meant to set.
      throw std::invalid_argument(location + "unknown parameter " + segs[0]);
    }
    const bool is_list = segs[0] == "element_set" || segs[0] == "solute_element_set" ||
        segs[0] == "solute_number_set" || segs[0] == "vacancy_trajectory";
    if (!is_list && segs.size() != 2) {
      throw std::invalid_argument(location + "parameter " + segs[0] + " takes one value but " +
                                  std::to_string(segs.size() - 1) + " were given");
    }
  }
  ifs.close();
}
}    // namespace api
