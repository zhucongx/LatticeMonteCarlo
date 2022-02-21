#include "ChainKmcOmp.h"

#include <utility>
#include <chrono>
#include <utility>
namespace kmc {
ChainKmcOmp::ChainKmcOmp(cfg::Config config,
                         unsigned long long int log_dump_steps,
                         unsigned long long int config_dump_steps,
                         unsigned long long int maximum_number,
                         const std::set<Element> &type_set,
                         unsigned long long int steps,
                         double energy,
                         double time,
                         const std::string &json_parameters_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      steps_(steps),
      energy_(energy),
      time_(time),
      vacancy_index_(cfg::GetVacancyAtomIndex(config_)),
      previous_j_(config_.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_)[0]),
      energy_predictor_(json_parameters_filename,
                        config_, type_set),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  event_list_.resize(kFirstEventListSize);
}
ChainKmcOmp::~ChainKmcOmp() = default;
void ChainKmcOmp::Dump(std::ofstream &ofs) {
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << time_ << '\t' << energy_ << '\t' << one_step_barrier_ << '\t'
        << one_step_energy_change_ << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteCfg(std::to_string(steps_) + ".cfg", true);
  }
}
// Get the energy change and the probability from j to k, pjk by the reference.
// And return the index of j and k. Only applied to 12 sub-primary processes.
// And then pass to others, now we will have 12 different numbers for 12 different second groups.

} // namespace kmc