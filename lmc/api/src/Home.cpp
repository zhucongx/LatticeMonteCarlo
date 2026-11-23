#include "Home.h"
#include <mpi.h>
#include <iomanip>

bool IsRoot() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
}

namespace api {
namespace {
// Stops before anything expensive is built when a method is started with a number of MPI processes it cannot use.
void RequireNumProcesses(const std::string &method, const int required_size) {
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != required_size) {
    throw std::invalid_argument("simulation_method " + method + " must be run with " + std::to_string(required_size) +
                                " MPI process(es), but " + std::to_string(world_size) + " started");
  }
}
}    // namespace

void Print(const Parameter &parameter) {
  if (!IsRoot()) {
    return;
  }
  for (const auto &warning: parameter.warnings_) {
    std::cerr << "Warning: " << warning << std::endl;
  }
  std::cout << "Parameters" << std::endl;
  std::cout << "simulation_method: " << parameter.method << std::endl;
  if (parameter.method == "KineticMcFirstMpi" || parameter.method == "KineticMcFirstSerial" ||
      parameter.method == "KineticMcSecondOmpi") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_ << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "map_filename: " << parameter.map_filename_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(),
              parameter.element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_ << std::endl;
    std::cout << "temperature: " << parameter.temperature_ << std::endl;
    std::cout << "restart_steps: " << parameter.restart_steps_ << std::endl;
    std::cout << "restart_energy: " << parameter.restart_energy_ << std::endl;
    std::cout << "restart_time: " << parameter.restart_time_ << std::endl;
    std::cout << "time_temperature_filename: " << parameter.time_temperature_filename_ << std::endl;
    std::cout << std::boolalpha;
    std::cout << "early_stop: " << parameter.early_stop_ << std::endl;
    std::cout << "solute_disp: " << parameter.solute_disp_ << std::endl;
    std::cout << std::noboolalpha;
    std::cout << "seed: " << parameter.seed_ << std::endl;
  } else if (parameter.method == "SimulatedAnnealing") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_ << std::endl;
    std::cout << "factor: " << parameter.factor_ << std::endl;
    std::cout << "solvent_element: " << parameter.solvent_element_ << std::endl;
    std::cout << "solute_element_set: ";
    std::copy(parameter.solute_element_set_.begin(),
              parameter.solute_element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "solute_number_set: ";
    std::transform(parameter.solute_number_set_.begin(),
                   parameter.solute_number_set_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto number) {
                     return std::to_string(number);
                   });
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "initial_temperature: " << parameter.initial_temperature_ << std::endl;
    std::cout << "seed: " << parameter.seed_ << std::endl;
  } else if (parameter.method == "CanonicalMcSerial" || parameter.method == "CanonicalMcOmp") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_ << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "map_filename: " << parameter.map_filename_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(),
              parameter.element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_ << std::endl;
    std::cout << "temperature: " << parameter.temperature_ << std::endl;
    std::cout << "restart_steps: " << parameter.restart_steps_ << std::endl;
    std::cout << "restart_energy: " << parameter.restart_energy_ << std::endl;
    std::cout << "seed: " << parameter.seed_ << std::endl;
  } else if (parameter.method == "Ansys" || parameter.method == "Reformat") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_ << std::endl;
    std::cout << "solvent_element: " << parameter.solvent_element_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(),
              parameter.element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "initial_steps: " << parameter.initial_steps_ << std::endl;
    std::cout << "increment_steps: " << parameter.increment_steps_ << std::endl;
    std::cout << "smallest_cluster_criteria: " << parameter.smallest_cluster_criteria_ << std::endl;
    std::cout << "solvent_bond_criteria: " << parameter.solvent_bond_criteria_ << std::endl;
    std::cout << "log_type: " << parameter.log_type_ << std::endl;
    std::cout << "config_type: " << parameter.config_type_ << std::endl;
  }
}

void Run(const Parameter &parameter) {
  // The two MPI variants evaluate one of the twelve jumps of the vacancy on each rank. Everything else runs one
  // trajectory, or one analysis, in a single process.
  if (parameter.method == "KineticMcFirstMpi" || parameter.method == "KineticMcSecondOmpi") {
    RequireNumProcesses(parameter.method, static_cast<int>(constants::kNumFirstNearestNeighbors));
  } else {
    RequireNumProcesses(parameter.method, 1);
  }
  if (parameter.method == "KineticMcFirstMpi") {
    auto kinetic_mc_first_mpi = api::BuildKineticMcFirstMpiFromParameter(parameter);
    kinetic_mc_first_mpi.Simulate();
  } else if (parameter.method == "KineticMcFirstSerial") {
    auto kinetic_mc_first_serial = api::BuildKineticMcFirstSerialFromParameter(parameter);
    kinetic_mc_first_serial.Simulate();
  } else if (parameter.method == "KineticMcSecondOmpi") {
    auto kinetic_mc_second_ompi = api::BuildKineticMcSecondOmpiFromParameter(parameter);
    kinetic_mc_second_ompi.Simulate();
  } else if (parameter.method == "Ansys") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunAnsys();
  } else if (parameter.method == "Reformat") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunReformat();
  } else if (parameter.method == "SimulatedAnnealing") {
    auto simulated_annealing = api::BuildSimulatedAnnealingFromParameter(parameter);
    simulated_annealing.Simulate();
  } else if (parameter.method == "CanonicalMcSerial") {
    auto canonical_mc_serial = api::BuildCanonicalMcSerialFromParameter(parameter);
    canonical_mc_serial.Simulate();
  } else if (parameter.method == "CanonicalMcOmp") {
    auto canonical_mc_omp = api::BuildCanonicalMcOmpFromParameter(parameter);
    canonical_mc_omp.Simulate();
  } else {
    throw std::invalid_argument("No such method: " + parameter.method);
  }
}

mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadConfig(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  if (IsRoot()) {
    std::cout << "Finish config reading. Start KMC." << std::endl;
  }
  return mc::KineticMcFirstMpi{config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.thermodynamic_averaging_steps_,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.restart_time_,
                               parameter.temperature_,
                               element_set,
                               parameter.json_coefficients_filename_,
                               parameter.time_temperature_filename_,
                               parameter.early_stop_,
                               parameter.solute_disp_,
                               parameter.seed_};
}

mc::KineticMcFirstSerial BuildKineticMcFirstSerialFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadConfig(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  if (IsRoot()) {
    std::cout << "Finish config reading. Start KMC." << std::endl;
  }
  return mc::KineticMcFirstSerial{config,
                                  parameter.log_dump_steps_,
                                  parameter.config_dump_steps_,
                                  parameter.maximum_steps_,
                                  parameter.thermodynamic_averaging_steps_,
                                  parameter.restart_steps_,
                                  parameter.restart_energy_,
                                  parameter.restart_time_,
                                  parameter.temperature_,
                                  element_set,
                                  parameter.json_coefficients_filename_,
                                  parameter.time_temperature_filename_,
                                  parameter.early_stop_,
                                  parameter.solute_disp_,
                                  parameter.seed_};
}

mc::KineticMcSecondOmpi BuildKineticMcSecondOmpiFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadConfig(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  if (IsRoot()) {
    std::cout << "Finish config reading. Start KMC." << std::endl;
  }
  return mc::KineticMcSecondOmpi{config,
                                 parameter.log_dump_steps_,
                                 parameter.config_dump_steps_,
                                 parameter.maximum_steps_,
                                 parameter.thermodynamic_averaging_steps_,
                                 parameter.restart_steps_,
                                 parameter.restart_energy_,
                                 parameter.restart_time_,
                                 parameter.temperature_,
                                 element_set,
                                 parameter.json_coefficients_filename_,
                                 parameter.time_temperature_filename_,
                                 parameter.early_stop_,
                                 parameter.solute_disp_,
                                 parameter.seed_};
}

mc::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter) {
  std::map<Element, size_t> solute_atom_count;
  for (size_t i = 0; i < parameter.solute_element_set_.size(); ++i) {
    solute_atom_count.insert(
        std::make_pair(Element(parameter.solute_element_set_[i]), parameter.solute_number_set_[i]));
  }

  return mc::SimulatedAnnealing{{parameter.factor_, parameter.factor_, parameter.factor_},
                                Element(parameter.solvent_element_),
                                solute_atom_count,
                                parameter.log_dump_steps_,
                                parameter.config_dump_steps_,
                                parameter.maximum_steps_,
                                parameter.initial_temperature_,
                                parameter.json_coefficients_filename_,
                                parameter.seed_};
}

mc::CanonicalMcSerial BuildCanonicalMcSerialFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadConfig(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  if (IsRoot()) {
    std::cout << "Finish config reading. Start CMC." << std::endl;
  }
  return mc::CanonicalMcSerial{config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.thermodynamic_averaging_steps_,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.temperature_,
                               element_set,
                               parameter.json_coefficients_filename_,
                               parameter.seed_};
}

mc::CanonicalMcOmp BuildCanonicalMcOmpFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadConfig(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  if (IsRoot()) {
    std::cout << "Finish config reading. Start CMC." << std::endl;
  }
  return mc::CanonicalMcOmp{config,
                            parameter.log_dump_steps_,
                            parameter.config_dump_steps_,
                            parameter.maximum_steps_,
                            parameter.thermodynamic_averaging_steps_,
                            parameter.restart_steps_,
                            parameter.restart_energy_,
                            parameter.temperature_,
                            element_set,
                            parameter.json_coefficients_filename_,
                            parameter.seed_};
}

ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  return ansys::Traverse{parameter.initial_steps_,
                         parameter.increment_steps_,
                         Element(parameter.solvent_element_),
                         element_set,
                         parameter.smallest_cluster_criteria_,
                         parameter.solvent_bond_criteria_,
                         parameter.escape_temperature_,
                         parameter.json_coefficients_filename_,
                         parameter.log_type_,
                         parameter.config_type_};
}
}    // namespace api
