#include "Home.h"
namespace api {
void Print(const Parameter &parameter) {
  std::cout << "Parameters" << std::endl;
  std::cout << "simulation_method: " << parameter.method << std::endl;
  if (parameter.method == "KineticMcFirstMpi" || parameter.method == "KineticMcChainOmp"
      || parameter.method == "KineticMcChainMpi" || parameter.method == "KineticMcChainOmpi") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "map_filename: " << parameter.map_filename_ << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "temperature: " << parameter.temperature_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(), parameter.element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "restart_steps: " << parameter.restart_steps_ << std::endl;
    std::cout << "restart_energy: " << parameter.restart_energy_ << std::endl;
    std::cout << "restart_time: " << parameter.restart_time_ << std::endl;
  } else if (parameter.method == "SimulatedAnnealing") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "factor: " << parameter.factor_ << std::endl;
    std::cout << "solvent_element: " << parameter.solvent_element_ << std::endl;
    std::cout << "solute_element_set: ";
    std::copy(parameter.solute_element_set_.begin(), parameter.solute_element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "solute_number_set: ";
    std::transform(parameter.solute_number_set_.begin(), parameter.solute_number_set_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto number) { return std::to_string(number); });
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "early_stop_steps: " << parameter.early_stop_steps_ << std::endl;
    std::cout << "initial_temperature: " << parameter.initial_temperature_ << std::endl;
  } else if (parameter.method == "CanonicalMc") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(), parameter.element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "temperature: " << parameter.temperature_ << std::endl;
  } else if (parameter.method == "CanonicalMcStepT"
      || parameter.method == "SemiGrandCanonicalMcStepT") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(), parameter.element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "initial_temperature: " << parameter.initial_temperature_ << std::endl;
    std::cout << "decrement_temperature: " << parameter.decrement_temperature_ << std::endl;

  } else if (parameter.method == "FindCluster" || parameter.method == "ShortRangeOrder"
      || parameter.method == "Reformat") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "solvent_element: " << parameter.solvent_element_ << std::endl;
    std::cout << "element_set: ";
    std::copy(parameter.element_set_.begin(), parameter.element_set_.end(),
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
  if (parameter.method == "KineticMcFirstMpi") {
    auto kinetic_mc_first_mpi = api::BuildKineticMcFirstMpiFromParameter(parameter);
    kinetic_mc_first_mpi.Simulate();
  } else if (parameter.method == "KineticMcChainOmp") {
    auto kinetic_mc_chain_omp = api::BuildKineticMcChainOmpFromParameter(parameter);
    kinetic_mc_chain_omp.Simulate();
  } else if (parameter.method == "KineticMcChainMpi") {
    auto kinetic_mc_chain_mpi = api::BuildKineticMcChainMpiFromParameter(parameter);
    kinetic_mc_chain_mpi.Simulate();
  } else if (parameter.method == "KineticMcChainOmpi") {
    auto kinetic_mc_chain_ompi = api::BuildKineticMcChainOmpiFromParameter(parameter);
    kinetic_mc_chain_ompi.Simulate();
  } else if (parameter.method == "Cluster") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunCluster();
  } else if (parameter.method == "ShortRangeOrder") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunShortRangeOrder();
  } else if (parameter.method == "Reformat") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunReformat();
  } else if (parameter.method == "SimulatedAnnealing") {
    auto simulated_annealing = api::BuildSimulatedAnnealingFromParameter(parameter);
    simulated_annealing.Simulate();
  } else if (parameter.method == "CanonicalMc") {
    auto canonical_mc = api::BuildCanonicalMcFromParameter(parameter);
    canonical_mc.Simulate();
  } else if (parameter.method == "CanonicalMcStepT") {
    auto canonical_mc_step_t = api::BuildCanonicalMcStepTFromParameter(parameter);
    canonical_mc_step_t.Simulate();
  } else if (parameter.method == "SemiGrandCanonicalMcStepT") {
    auto semi_grand_canonical_mc_step_t =
        api::BuildSemiGrandCanonicalMcStepTFromParameter(parameter);
    semi_grand_canonical_mc_step_t.Simulate();
  } else {
    std::cout << "No such method: " << parameter.method << std::endl;
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
  std::cout << "Finish config reading. Start KMC." << std::endl;
  return mc::KineticMcFirstMpi{config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.temperature_,
                               element_set,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.restart_time_,
                               parameter.json_coefficients_filename_};
}
mc::KineticMcChainMpi BuildKineticMcChainMpiFromParameter(const Parameter &parameter) {
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
  std::cout << "Finish config reading. Start KMC." << std::endl;
  return mc::KineticMcChainMpi{config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.temperature_,
                               element_set,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.restart_time_,
                               parameter.json_coefficients_filename_};
}
mc::KineticMcChainOmp BuildKineticMcChainOmpFromParameter(const Parameter &parameter) {
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
  std::cout << "Finish config reading. Start kMC." << std::endl;
  return mc::KineticMcChainOmp{config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.temperature_,
                               element_set,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.restart_time_,
                               parameter.json_coefficients_filename_};
}
mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter &parameter) {
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
  std::cout << "Finish config reading. Start KMC." << std::endl;
  return mc::KineticMcChainOmpi{config,
                                parameter.log_dump_steps_,
                                parameter.config_dump_steps_,
                                parameter.maximum_steps_,
                                parameter.temperature_,
                                element_set,
                                parameter.restart_steps_,
                                parameter.restart_energy_,
                                parameter.restart_time_,
                                parameter.json_coefficients_filename_};
}
ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter) {
  std::map<Element, size_t> solute_atom_count;
  for (size_t i = 0; i < parameter.solute_element_set_.size(); ++i) {
    solute_atom_count.insert(std::make_pair(Element(parameter.solute_element_set_[i]),
                                            parameter.solute_number_set_[i]));
  }

  return ansys::SimulatedAnnealing{
      {parameter.factor_, parameter.factor_, parameter.factor_},
      Element(parameter.solvent_element_),
      solute_atom_count,
      parameter.log_dump_steps_,
      parameter.config_dump_steps_,
      parameter.maximum_steps_,
      parameter.early_stop_steps_,
      parameter.initial_temperature_,
      parameter.json_coefficients_filename_};

}
mc::CanonicalMc BuildCanonicalMcFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  auto config = cfg::Config::ReadConfig(parameter.config_filename_);
  std::cout << "Finish config reading. Start CMC." << std::endl;
  return mc::CanonicalMc{
      config,
      element_set,
      parameter.log_dump_steps_,
      parameter.config_dump_steps_,
      parameter.maximum_steps_,
      parameter.temperature_,
      parameter.json_coefficients_filename_};
}
mc::CanonicalMcStepT BuildCanonicalMcStepTFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  auto config = cfg::Config::ReadConfig(parameter.config_filename_);
  std::cout << "Finish config reading. Start CMC." << std::endl;
  return mc::CanonicalMcStepT{
      config,
      element_set,
      parameter.log_dump_steps_,
      parameter.config_dump_steps_,
      parameter.maximum_steps_,
      parameter.initial_temperature_,
      parameter.decrement_temperature_,
      parameter.json_coefficients_filename_};
}
mc::SemiGrandCanonicalMcStepT BuildSemiGrandCanonicalMcStepTFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  auto config = cfg::Config::ReadConfig(parameter.config_filename_);
  std::cout << "Finish config reading. Start SGCMC." << std::endl;
  return mc::SemiGrandCanonicalMcStepT{
      config,
      element_set,
      parameter.log_dump_steps_,
      parameter.config_dump_steps_,
      parameter.maximum_steps_,
      parameter.initial_temperature_,
      parameter.decrement_temperature_,
      parameter.json_coefficients_filename_};
}
ansys::Iterator BuildIteratorFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  return ansys::Iterator{parameter.initial_steps_, parameter.increment_steps_,
                         Element(parameter.solvent_element_), element_set,
                         parameter.smallest_cluster_criteria_,
                         parameter.solvent_bond_criteria_,
                         parameter.json_coefficients_filename_,
                         parameter.log_type_,
                         parameter.config_type_};
}
} // api