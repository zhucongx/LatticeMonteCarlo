#include "ChainKmcOmp.h"

#include <utility>
#include <chrono>
#include <utility>
namespace kmc {

ChainKmcOmp::ChainKmcOmp(cfg::Config config,
                         unsigned long long int log_dump_steps,
                         unsigned long long int config_dump_steps,
                         unsigned long long int maximum_number,
                         double temperature,
                         const std::set<Element> &type_set,
                         unsigned long long int steps,
                         double energy,
                         double time,
                         const std::string &json_parameters_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      coefficient_(1 / kBoltzmannConstant / temperature),
      steps_(steps),
      energy_(energy),
      time_(time),
      vacancy_index_(cfg::GetVacancyAtomIndex(config_)),
      previous_j_(config_.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_)[0]),
      energy_predictor_(json_parameters_filename,
                        config_, type_set),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using upto" << omp_get_num_threads() << " threads." << std::endl;
  }
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
void ChainKmcOmp::Simulate() {

}
// Get the energy change and the probability from j to k, pjk by the reference.
// And return the index of j and k. Only applied to 12 threads.
void ChainKmcOmp::BuildEventList() {
  total_rate_k_ = 0;
  const auto i_indexes = config_.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_);
#pragma omp parallel default(none) shared(i_indexes) reduction(+: total_rate_k_)
  {
#pragma omp for
    for (size_t it = 0; it < kFirstEventListSize; ++it) {
      auto i_index = i_indexes[it];
      auto event_i = JumpEvent(
          {vacancy_index_, i_index},
          energy_predictor_.GetBarrierAndDiffFromAtomIdPair(config_,
                                                            {vacancy_index_,
                                                             i_index}),
          coefficient_);
      total_rate_k_ += event_i.GetForwardRate();
      event_list_[it] = std::move(event_i);
    }
  }
}
double ChainKmcOmp::BuildEventListParallel() {
  BuildEventList();
#pragma omp parallel default(none)
  {
#pragma omp for collapse(2)
    for (size_t it1 = 0; it1 < kFirstEventListSize; ++it1) {
      for (size_t it2 = 0; it2 < kSecondEventListSize; ++it2) {

      }
    }
  }
  return 0;
}
} // namespace kmc