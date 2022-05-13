#include "Home.h"
namespace api {
kmc::FirstKmcMpi BuildFirstKmcMpiFromParameter(const Parameter &parameter) {
  std::set<Element> ele_set;
  for (const auto &element_string: parameter.element_symbols_) {
    ele_set.insert(Element(element_string));
  }
  return kmc::FirstKmcMpi{cfg::Config::ReadCfg(parameter.config_filename_),
                          parameter.log_dump_steps_,
                          parameter.config_dump_steps_,
                          parameter.maximum_number_,
                          parameter.temperature_,
                          ele_set,
                          parameter.restart_steps_,
                          parameter.restart_energy_,
                          parameter.restart_time_,
                          parameter.json_coefficients_filename_};
}
kmc::ChainKmcMpi BuildChainKmcMpiFromParameter(const Parameter &parameter) {
  std::set<Element> ele_set;
  for (const auto &element_string: parameter.element_symbols_) {
    ele_set.insert(Element(element_string));
  }
  return kmc::ChainKmcMpi{cfg::Config::ReadCfg(parameter.config_filename_),
                          parameter.log_dump_steps_,
                          parameter.config_dump_steps_,
                          parameter.maximum_number_,
                          parameter.temperature_,
                          ele_set,
                          parameter.restart_steps_,
                          parameter.restart_energy_,
                          parameter.restart_time_,
                          parameter.json_coefficients_filename_};
}

} // namespace api