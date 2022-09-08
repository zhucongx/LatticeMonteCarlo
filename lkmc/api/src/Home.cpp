#include "Home.h"
namespace api {
kmc::FirstKmcMpi BuildFirstKmcMpiFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadCfg(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  std::cout << "Finish config reading. Start kMC." << std::endl;
  return kmc::FirstKmcMpi{config,
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
kmc::ChainKmcMpi BuildChainKmcMpiFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }
  cfg::Config config;
  if (parameter.map_filename_.empty()) {
    config = cfg::Config::ReadCfg(parameter.config_filename_);
    config.ReassignLatticeVector();
  } else {
    config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }
  std::cout << "Finish config reading. Start kMC." << std::endl;
  return kmc::ChainKmcMpi{config,
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
      parameter.temperature_,
      parameter.json_coefficients_filename_};

}
ansys::CanonicalMC BuildCanonicalMCFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  auto config = cfg::Config::ReadCfg(parameter.config_filename_);
  return ansys::CanonicalMC{
      config,
      element_set,
      parameter.log_dump_steps_,
      parameter.config_dump_steps_,
      parameter.maximum_steps_,
      parameter.temperature_,
      parameter.json_coefficients_filename_};

}
ansys::Iterator BuildIteratorFromParameter(const Parameter &parameter) {
  std::set<Element> element_set;
  for (const auto &element_string: parameter.element_set_) {
    element_set.insert(Element(element_string));
  }

  return ansys::Iterator{parameter.initial_steps_,
                         parameter.increment_steps_,
                         Element(parameter.solvent_element_),
                         element_set,
                         parameter.smallest_cluster_criteria_,
                         parameter.json_coefficients_filename_};
}
} // namespace api