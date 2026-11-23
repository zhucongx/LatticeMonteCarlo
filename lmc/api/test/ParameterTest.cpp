#include "Parameter.h"
#include "cfg/test/TemporaryPath.hpp"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace api {
namespace {
const std::string kKineticMc = "simulation_method KineticMcFirstSerial\n"
                               "json_coefficients_filename quartic_coefficients.json\n"
                               "config_filename start.cfg\n"
                               "log_dump_steps 100\n"
                               "config_dump_steps 1000\n"
                               "maximum_steps 5000\n"
                               "temperature 450\n"
                               "element_set Al Mg Zn\n";

const std::string kSimulatedAnnealing = "simulation_method SimulatedAnnealing\n"
                                        "json_coefficients_filename quartic_coefficients.json\n"
                                        "factor 6\n"
                                        "solvent_element Al\n"
                                        "solute_element_set Mg Zn\n"
                                        "solute_number_set 10 12\n"
                                        "log_dump_steps 100\n"
                                        "maximum_steps 5000\n"
                                        "initial_temperature 700\n";

const std::string kAnalysis = "simulation_method Ansys\n"
                              "log_type kinetic_mc\n"
                              "config_type config\n"
                              "json_coefficients_filename quartic_coefficients.json\n"
                              "solvent_element Al\n"
                              "element_set Al Mg Zn\n"
                              "initial_steps 0\n"
                              "increment_steps 10000\n"
                              "smallest_cluster_criteria 4\n"
                              "solvent_bond_criteria 3\n"
                              "escape_temperature 300\n";

// Writes the text to a parameter file of its own and reads it.
Parameter Read(const std::string &content) {
  const test::TemporaryFile file(content, ".txt");
  return Parameter(file.GetName());
}

std::string Replace(std::string text, const std::string &from, const std::string &to) {
  const auto position = text.find(from);
  EXPECT_NE(position, std::string::npos) << "'" << from << "' is not part of the parameter text";
  if (position != std::string::npos) {
    text.replace(position, from.size(), to);
  }
  return text;
}

// The text without the line that sets `key`.
std::string Without(const std::string &text, const std::string &key) {
  const auto begin = text.find(key + " ");
  EXPECT_NE(begin, std::string::npos) << key;
  return text.substr(0, begin) + text.substr(text.find('\n', begin) + 1);
}
}    // namespace

TEST(Parameter, ReadsAKineticMcFile) {
  const auto parameter = Read(kKineticMc);
  EXPECT_EQ(parameter.method, "KineticMcFirstSerial");
  EXPECT_EQ(parameter.json_coefficients_filename_, "quartic_coefficients.json");
  EXPECT_EQ(parameter.config_filename_, "start.cfg");
  EXPECT_EQ(parameter.log_dump_steps_, 100u);
  EXPECT_EQ(parameter.config_dump_steps_, 1000u);
  EXPECT_EQ(parameter.maximum_steps_, 5000u);
  EXPECT_DOUBLE_EQ(parameter.temperature_, 450.0);
  EXPECT_EQ(parameter.element_set_, (std::vector<std::string>{"Al", "Mg", "Zn"}));
  EXPECT_FALSE(parameter.early_stop_);
  EXPECT_TRUE(parameter.warnings_.empty());
}

// A line of blanks used to crash, a tab used to hide the value and a carriage return used to end up in it.
TEST(Parameter, AcceptsTabsBlankLinesCommentsAndWindowsLineEndings) {
  std::string content = "# comment\n\n \t \n" + Replace(kKineticMc, "temperature 450", "  temperature\t450   # Kelvin");
  for (auto position = content.find('\n'); position != std::string::npos; position = content.find('\n', position + 2)) {
    content.replace(position, 1, "\r\n");
  }
  const auto parameter = Read(content);
  EXPECT_EQ(parameter.method, "KineticMcFirstSerial");
  EXPECT_DOUBLE_EQ(parameter.temperature_, 450.0);
  EXPECT_EQ(parameter.element_set_, (std::vector<std::string>{"Al", "Mg", "Zn"}));
}

// A misspelled key used to leave its parameter at zero without any message.
TEST(Parameter, RejectsUnknownKey) {
  EXPECT_THROW(Read(Replace(kKineticMc, "temperature 450", "temperture 450")), std::invalid_argument);
  EXPECT_THROW(Read(kKineticMc + "no_such_parameter 1\n"), std::invalid_argument);
}

TEST(Parameter, RejectsKeyWithoutValueOrWithTooManyValues) {
  EXPECT_THROW(Read(Replace(kKineticMc, "temperature 450", "temperature")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "temperature 450", "temperature 450 K")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "start.cfg", "my start.cfg")), std::invalid_argument);
}

// std::stoull and std::stod accept a prefix: "1e3" was 1 step, "-5" a huge number and "450K" was 450.
TEST(Parameter, RejectsMalformedNumbers) {
  EXPECT_THROW(Read(Replace(kKineticMc, "maximum_steps 5000", "maximum_steps 5e3")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "maximum_steps 5000", "maximum_steps -5000")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "maximum_steps 5000", "maximum_steps 5000.0")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "temperature 450", "temperature 450K")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "temperature 450", "temperature nan")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kSimulatedAnnealing, "10 12", "10 twelve")), std::invalid_argument);
  EXPECT_DOUBLE_EQ(Read(Replace(kKineticMc, "temperature 450", "temperature 4.5e2")).temperature_, 450.0);
}

// Everything except the exact word "true" used to be read as false.
TEST(Parameter, ReadsBooleans) {
  const auto parameter = Read(kKineticMc + "early_stop 0\nsolute_disp 1\n");
  EXPECT_FALSE(parameter.early_stop_);
  EXPECT_TRUE(parameter.solute_disp_);
  EXPECT_TRUE(Read(kKineticMc + "early_stop True\n").early_stop_);
  EXPECT_THROW(Read(kKineticMc + "early_stop yes\n"), std::invalid_argument);
}

TEST(Parameter, RequiresTheParametersOfTheMethod) {
  for (const std::string key: {"simulation_method", "json_coefficients_filename", "config_filename", "log_dump_steps",
                               "config_dump_steps", "maximum_steps", "temperature", "element_set"}) {
    EXPECT_THROW(Read(Without(kKineticMc, key)), std::invalid_argument) << key;
  }
  for (const std::string key: {"factor", "solvent_element", "solute_element_set", "solute_number_set",
                               "log_dump_steps", "maximum_steps", "initial_temperature"}) {
    EXPECT_THROW(Read(Without(kSimulatedAnnealing, key)), std::invalid_argument) << key;
  }
  for (const std::string key: {"log_type", "config_type", "solvent_element", "element_set", "increment_steps",
                               "smallest_cluster_criteria", "solvent_bond_criteria", "escape_temperature"}) {
    EXPECT_THROW(Read(Without(kAnalysis, key)), std::invalid_argument) << key;
  }
  EXPECT_NO_THROW(Read(kSimulatedAnnealing));
  EXPECT_NO_THROW(Read(kAnalysis));
}

TEST(Parameter, AcceptsMapFileInsteadOfConfigFile) {
  const auto parameter = Read(Replace(kKineticMc, "config_filename start.cfg", "map_filename map0.txt"));
  EXPECT_EQ(parameter.map_filename_, "map0.txt");
}

// With a time-temperature table the temperature comes from the table.
TEST(Parameter, TemperatureIsOptionalWithATimeTemperatureTable) {
  EXPECT_NO_THROW(Read(Without(kKineticMc, "temperature") + "time_temperature_filename time_temperature.dat\n"));
}

TEST(Parameter, RejectsValuesTheRunCannotUse) {
  EXPECT_THROW(Read(Replace(kKineticMc, "temperature 450", "temperature 0")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "log_dump_steps 100", "log_dump_steps 0")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "config_dump_steps 1000", "config_dump_steps 0")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kAnalysis, "increment_steps 10000", "increment_steps 0")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kAnalysis, "log_type kinetic_mc", "log_type kmc")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kSimulatedAnnealing, "solute_number_set 10 12", "solute_number_set 10")),
               std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "KineticMcFirstSerial", "KineticMcFirst")), std::invalid_argument);
  // The names KineticMcFirstSerial and KineticMcSecondOmpi had before are not accepted any more.
  EXPECT_THROW(Read(Replace(kKineticMc, "KineticMcFirstSerial", "KineticMcFirstOmp")), std::invalid_argument);
  EXPECT_THROW(Read(Replace(kKineticMc, "KineticMcFirstSerial", "KineticMcChainOmpi")), std::invalid_argument);
}

TEST(Parameter, AcceptsTheThreeKineticMcMethods) {
  for (const std::string method: {"KineticMcFirstSerial", "KineticMcFirstMpi", "KineticMcSecondOmpi"}) {
    EXPECT_EQ(Read(Replace(kKineticMc, "KineticMcFirstSerial", method)).method, method);
  }
}

TEST(Parameter, WarnsAboutObsoleteAndRepeatedKeys) {
  const auto obsolete = Read(kKineticMc + "vacancy_trajectory 0 0 0\n");
  ASSERT_EQ(obsolete.warnings_.size(), 1u);
  EXPECT_NE(obsolete.warnings_[0].find("vacancy_trajectory"), std::string::npos);

  const auto repeated = Read(kKineticMc + "temperature 500\n");
  ASSERT_EQ(repeated.warnings_.size(), 1u);
  EXPECT_NE(repeated.warnings_[0].find("temperature"), std::string::npos);
  EXPECT_DOUBLE_EQ(repeated.temperature_, 500.0);
}

// Earlier versions could rescale the time of a KMC run. A file that switched this off still loads, a file that
// switched it on asks for other times than the program writes now.
TEST(Parameter, RateCorrectorOfEarlierVersionsCanOnlyBeOff) {
  const auto parameter = Read(kKineticMc + "rate_corrector false\n");
  ASSERT_EQ(parameter.warnings_.size(), 1u);
  EXPECT_NE(parameter.warnings_[0].find("rate_corrector"), std::string::npos);
  EXPECT_THROW(Read(kKineticMc + "rate_corrector true\n"), std::invalid_argument);
}

TEST(Parameter, SeedIsTakenFromTheFile) {
  EXPECT_EQ(Read(kKineticMc + "seed 12345\n").seed_, 12345u);
  EXPECT_THROW(Read(kKineticMc + "seed -1\n"), std::invalid_argument);
}
}    // namespace api
