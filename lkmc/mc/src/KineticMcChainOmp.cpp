#include "KineticMcChainOmp.h"

#include <utility>
#include <chrono>
#include <utility>
namespace mc {
//  j -> k -> i ->l
//       |
// current position
// Todo: calculate ik event rather than kj event
KineticMcChainOmp::KineticMcChainOmp(const cfg::Config &config,
                                     unsigned long long int log_dump_steps,
                                     unsigned long long int config_dump_steps,
                                     unsigned long long int maximum_steps,
                                     double temperature,
                                     const std::set<Element> &element_set,
                                     unsigned long long int restart_steps,
                                     double restart_energy,
                                     double restart_time,
                                     const std::string &json_coefficients_filename)
    : log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_steps_(maximum_steps),
      beta_(1.0 / constants::kBoltzmann / temperature),
      steps_(restart_steps),
      energy_(restart_energy),
      time_(restart_time),
      vacancy_index_(config.GetVacancyAtomIndex()),
      previous_j_(config.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_)[0]),
      energy_predictor_(json_coefficients_filename,
                        config, element_set, 100000),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  std::array<cfg::Config, kSecondEventListSize> config_sublist;
  config_list_.fill(config);
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
  }
}
KineticMcChainOmp::~KineticMcChainOmp() = default;
void KineticMcChainOmp::Dump(std::ofstream &ofs) const {
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << time_ << '\t' << energy_ << '\t' << one_step_barrier_ << '\t'
        << one_step_energy_change_ << '\t' << migrating_element_.GetString() << std::endl;
  }
  if (steps_ == 0) {
    config_list_[0].WriteLattice("lattice.txt");
    config_list_[0].WriteElement("element.txt");
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_list_[0].WriteMap("map" + std::to_string(steps_) + ".txt");
  }
}

// Get the energy change and the probability from j to k, pjk by the reference.
// Update event list and total rate of k and initial total i list. And find the indexes for l.
// Applied to 12 threads.
void KineticMcChainOmp::BuildEventKIListAndLIndexList() {
  total_rate_k_ = 0;
  const auto &config = config_list_[0];
  const auto &i_indexes = config.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_);

#pragma omp parallel for default(none) shared(i_indexes, config) reduction(+: total_rate_k_)
  for (size_t it = 0; it < kFirstEventListSize; ++it) {
    const auto i_index = i_indexes[it];
    auto event_k_i = JumpEvent(
        {vacancy_index_, i_index},
        energy_predictor_.GetBarrierAndDiffFromAtomIdPair(
            config, {vacancy_index_, i_index}),
        beta_);
    total_rate_k_ += event_k_i.GetForwardRate();
    // initial total rate i list
    total_rate_i_list_.at(it) = event_k_i.GetBackwardRate();
    // initial event list
    event_k_i_list_.at(it) = std::move(event_k_i);
    // update l_index_list_
    size_t ii = 0;
    for (const auto l_index: config.GetFirstNeighborsAtomIdVectorOfAtom(i_index)) {
      if (l_index == vacancy_index_) {
        continue;
      }
      l_index_list_[it * kSecondEventListSize + ii] = l_index;
      ++ii;
    }
  }
  for (auto &event_i: event_k_i_list_) {
    event_i.CalculateProbability(total_rate_k_);
  }
}

void KineticMcChainOmp::BuildTotalRateIList() {
#pragma omp parallel for default(none)
  for (size_t it = 0; it < kFirstEventListSize * kSecondEventListSize; ++it) {
    size_t it1 = it / kSecondEventListSize;

    const auto event_k_i = event_k_i_list_.at(it1);
    auto &config = config_list_.at(it);
    config.AtomJump(event_k_i.GetIdJumpPair());
    const auto l_index = l_index_list_.at(it);
    JumpEvent event_i_l({vacancy_index_, l_index},
                        energy_predictor_.GetBarrierAndDiffFromAtomIdPair(
                            config, {vacancy_index_, l_index}),
                        beta_);
    config.AtomJump(event_k_i.GetIdJumpPair());
    // get sum r_{k to l}
    auto r_i_l = event_i_l.GetForwardRate();
#pragma omp critical
    {
      total_rate_i_list_.at(it1) += r_i_l;
    }
  }
}

double KineticMcChainOmp::CalculateTime() {
  BuildEventKIListAndLIndexList();
  BuildTotalRateIList();
  double beta_bar_k = 0.0, beta_k = 0.0, gamma_bar_k_j = 0.0, gamma_k_j = 0.0,
      beta_k_j = 0.0, alpha_k_j = 0.0;
  std::array<double, kFirstEventListSize> beta_bar_k_i_list{}, beta_k_i_list{};
// #pragma omp parallel for default(none) reduction(+: beta_bar_k, beta_k, gamma_bar_k_j, gamma_k_j, beta_k_j, alpha_k_j)
  for (size_t it = 0; it < kFirstEventListSize; ++it) {
    const auto &event_k_i = event_k_i_list_.at(it);
    const auto probability_k_i = event_k_i.GetProbability();
    const auto probability_i_k = event_k_i.GetBackwardRate() / total_rate_i_list_.at(it);

    const auto beta_bar_k_i = probability_k_i * probability_i_k;
    const auto beta_k_i = probability_k_i * (1 - probability_i_k);
    beta_bar_k_i_list[it] = beta_bar_k_i;
    beta_k_i_list[it] = beta_k_i;
    beta_bar_k += beta_bar_k_i;
    beta_k += beta_k_i;

    double gamma_bar_k_j_helper = 0.0, gamma_k_j_helper = 0.0,
        beta_k_j_helper = 0.0, alpha_k_j_helper = 0.0;
    if (event_k_i.GetIdJumpPair().second == previous_j_) {
      beta_k_j_helper = beta_k_i;
      alpha_k_j_helper = probability_k_i;
    } else {
      gamma_bar_k_j_helper = beta_bar_k_i;
      gamma_k_j_helper = beta_k_i;
    }
    gamma_bar_k_j += gamma_bar_k_j_helper;
    gamma_k_j += gamma_k_j_helper;
    beta_k_j += beta_k_j_helper;
    alpha_k_j += alpha_k_j_helper;
  }
  const double one_over_one_minus_a_j = 1 / (1 - alpha_k_j);

  const double
      indirect_probability_k_j = one_over_one_minus_a_j * (gamma_bar_k_j / beta_k) * beta_k_j;
  for (size_t it = 0; it < kFirstEventListSize; ++it) {
    auto &event_k_i = event_k_i_list_.at(it);
    const auto beta_k_i = beta_k_i_list[it];
    const double
        indirect_probability_k_i = one_over_one_minus_a_j * (1 + gamma_bar_k_j / beta_k) * beta_k_i;
    if (event_k_i.GetIdJumpPair().second == previous_j_) {
      event_k_i.SetProbability(indirect_probability_k_j);
    } else {
      event_k_i.SetProbability(indirect_probability_k_i);
    }
  }
  // calculate relative and cumulative probability
  double cumulative_probability = 0.0;
  for (auto &event: event_k_i_list_) {
    cumulative_probability += event.GetProbability();
    event.SetCumulativeProbability(cumulative_probability);
  }
  double t = 1 / total_rate_k_ / constants::kPrefactor;
  double ts_numerator = 0.0, ts_j_numerator = 0.0;
  for (size_t it = 0; it < kFirstEventListSize; ++it) {
    const auto &event_k_i = event_k_i_list_.at(it);
    const auto total_rate_i = total_rate_i_list_.at(it);
    double t_i = 1 / total_rate_i / constants::kPrefactor;
    double ts_numerator_helper = (t + t_i) * beta_bar_k_i_list[it];
    ts_numerator += ts_numerator_helper;

    if (event_k_i.GetIdJumpPair().second == previous_j_) {
      ts_numerator_helper = 0;
    }
    ts_j_numerator += ts_numerator_helper;
  }
  double ts = ts_numerator / beta_bar_k;
  double ts_j = ts_j_numerator / gamma_bar_k_j;
  double t_2 = one_over_one_minus_a_j
      * (gamma_k_j * t + gamma_bar_k_j * (ts_j + t + beta_bar_k / beta_k * ts));

  return t_2;
}
size_t KineticMcChainOmp::SelectEvent() const {
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double random_number = distribution(generator_);
  auto it = std::lower_bound(event_k_i_list_.begin(),
                             event_k_i_list_.end(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.GetCumulativeProvability() < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == event_k_i_list_.cend()) {
    it--;
  }
  return static_cast<size_t>(std::distance(event_k_i_list_.begin(), it));
}

void KineticMcChainOmp::Simulate() {
  std::ofstream ofs("lkmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs << "steps\ttime\tenergy\tEa\tdE\ttype\n";
  ofs.precision(8);

  while (steps_ <= maximum_steps_) {
    Dump(ofs);

    one_step_time_change_ = CalculateTime();
    const JumpEvent &selected_event = event_k_i_list_[SelectEvent()];
    atom_id_jump_pair_ = selected_event.GetIdJumpPair();
    // update time and energy
    time_ += one_step_time_change_;
    one_step_energy_change_ = selected_event.GetEnergyChange();
    energy_ += one_step_energy_change_;
    one_step_barrier_ = selected_event.GetForwardBarrier();
    migrating_element_ = config_list_[0].GetAtomVector().at(atom_id_jump_pair_.second).GetElement();
#pragma omp parallel for default(none) shared(selected_event)
    for (auto &config: config_list_) {
      config.AtomJump(atom_id_jump_pair_);
    }
    previous_j_ = atom_id_jump_pair_.second;
    ++steps_;
  }
}
} // mc