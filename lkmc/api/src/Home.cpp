#include "Home.h"
namespace api {
kmc::FirstKmcMpi BuildFirstKmcMpiFromParameter(const Parameter &parameter) {
  std::set<Element> ele_set;
  for (const auto &element_string: parameter.element_set_) {
    ele_set.insert(Element(element_string));
  }
  return kmc::FirstKmcMpi{cfg::Config::ReadCfg(parameter.config_filename_),
                          parameter.log_dump_steps_,
                          parameter.config_dump_steps_,
                          parameter.maximum_steps_,
                          parameter.temperature_,
                          ele_set,
                          parameter.restart_steps_,
                          parameter.restart_energy_,
                          parameter.restart_time_,
                          parameter.json_coefficients_filename_};
}
kmc::ChainKmcMpi BuildChainKmcMpiFromParameter(const Parameter &parameter) {
  std::set<Element> ele_set;
  for (const auto &element_string: parameter.element_set_) {
    ele_set.insert(Element(element_string));
  }
  return kmc::ChainKmcMpi{cfg::Config::ReadCfg(parameter.config_filename_),
                          parameter.log_dump_steps_,
                          parameter.config_dump_steps_,
                          parameter.maximum_steps_,
                          parameter.temperature_,
                          ele_set,
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
} // namespace api