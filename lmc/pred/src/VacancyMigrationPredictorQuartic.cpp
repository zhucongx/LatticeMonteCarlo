#include "VacancyMigrationPredictorQuartic.h"
#include <Eigen/Dense>
#include <utility>
#include <boost/range/combine.hpp>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "Constants.hpp"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace pred {
namespace {
const std::vector<double>
    kDeClusterCounter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
}
VacancyMigrationPredictorQuartic::VacancyMigrationPredictorQuartic(const std::string &predictor_filename,
                                                                   const cfg::Config &reference_config,
                                                                   std::set<Element> element_set)
    : element_set_(std::move(element_set)),
      encoder_mmm_(element_set_, GetAverageClusterParametersMappingMMM(reference_config)),
      encoder_mm2_(element_set_, GetAverageClusterParametersMappingMM2(reference_config)),
      mapping_state_(FlattenClusterMapping(GetClusterParametersMappingStatePair(reference_config))),
      num_pair_slots_(reference_config.GetNumAtoms() * constants::kNumFirstNearestNeighbors),
      neighbor_slot_lookup_(reference_config.GetNumAtoms()) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);
  ordered_map_ = std::map<cfg::ElementCluster, int>(initialized_cluster_hashmap_.begin(),
                                                    initialized_cluster_hashmap_.end());
  std::vector<double> cluster_total_bonds;
  cluster_total_bonds.reserve(ordered_map_.size());
  for (const auto &entry : ordered_map_) {
    cluster_total_bonds.push_back(kDeClusterCounter.at(static_cast<size_t>(entry.first.GetLabel())));
  }
  state_cluster_indexer_ = ClusterIndexer(ordered_map_, std::move(cluster_total_bonds));

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  if (!ifs) {
    throw std::runtime_error("Cannot open " + predictor_filename);
  }
  nlohmann::json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters] : all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = JsonToEigenVector(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersQuartic{
        JsonToEigenVector(parameters.at("mu_x_mmm")),
        JsonToEigenVector(parameters.at("mu_x_mm2")),
        JsonToEigenVector(parameters.at("sigma_x_mmm")),
        JsonToEigenVector(parameters.at("sigma_x_mm2")),
        JsonToEigenMatrix(parameters.at("U_mmm")),
        JsonToEigenMatrix(parameters.at("U_mm2")),
        JsonToEigenVector(parameters.at("theta_D")),
        JsonToEigenVector(parameters.at("theta_Ks")),
        parameters.at("mu_D"),
        parameters.at("mu_Ks"),
        parameters.at("sigma_D"),
        parameters.at("sigma_Ks")
    };
  }
  CheckModelSize(predictor_filename, "Base.theta", static_cast<size_t>(base_theta_.size()),
                 state_cluster_indexer_.Size());
  const size_t encode_size_mmm = encoder_mmm_.Size();
  const size_t encode_size_mm2 = encoder_mm2_.Size();
  for (const auto &element : element_set_) {
    const auto it = element_parameters_hashmap_.find(element);
    if (it == element_parameters_hashmap_.end()) {
      throw std::runtime_error(predictor_filename + ": no coefficients for element " + element.GetString() +
                               " of element_set");
    }
    const auto &name = element.GetString();
    auto &parameters = it->second;
    CheckModelSize(predictor_filename, name + ".mu_x_mmm", static_cast<size_t>(parameters.mu_x_mmm.size()),
                   encode_size_mmm);
    CheckModelSize(predictor_filename, name + ".sigma_x_mmm", static_cast<size_t>(parameters.sigma_x_mmm.size()),
                   encode_size_mmm);
    CheckModelSize(predictor_filename, name + ".U_mmm (columns)", static_cast<size_t>(parameters.U_mmm.cols()),
                   encode_size_mmm);
    CheckModelSize(predictor_filename, name + ".theta_D", static_cast<size_t>(parameters.theta_D.size()),
                   static_cast<size_t>(parameters.U_mmm.rows()));
    CheckModelSize(predictor_filename, name + ".mu_x_mm2", static_cast<size_t>(parameters.mu_x_mm2.size()),
                   encode_size_mm2);
    CheckModelSize(predictor_filename, name + ".sigma_x_mm2", static_cast<size_t>(parameters.sigma_x_mm2.size()),
                   encode_size_mm2);
    CheckModelSize(predictor_filename, name + ".U_mm2 (columns)", static_cast<size_t>(parameters.U_mm2.cols()),
                   encode_size_mm2);
    CheckModelSize(predictor_filename, name + ".theta_Ks", static_cast<size_t>(parameters.theta_Ks.size()),
                   static_cast<size_t>(parameters.U_mm2.rows()));
    // The sizes are known to fit now. Fold the projection U, the coefficients theta, the scale of the encoding and
    // the averaging over the clusters of a group into one weight per cluster kind, so that an evaluation is a sum
    // over the clusters instead of an encoding and a product with a matrix of 150 rows.
    const Eigen::VectorXd weight_D =
        (parameters.U_mmm.transpose() * parameters.theta_D).cwiseQuotient(parameters.sigma_x_mmm);
    const Eigen::VectorXd weight_Ks =
        (parameters.U_mm2.transpose() * parameters.theta_Ks).cwiseQuotient(parameters.sigma_x_mm2);
    parameters.cluster_weight_D = encoder_mmm_.ScaleWeights(weight_D);
    parameters.cluster_weight_Ks = encoder_mm2_.ScaleWeights(weight_Ks);
    parameters.weight_dot_mu_D = weight_D.dot(parameters.mu_x_mmm);
    parameters.weight_dot_mu_Ks = weight_Ks.dot(parameters.mu_x_mm2);
  }
  site_bond_cluster_state_flat_.resize(num_pair_slots_);
  // The state, mmm and mm2 lists of a pair hold the same sites, sorted by three different rules in the frame of the
  // jump. If the sites of every pair lie in its frame where the sites of the first pair lie in its own, an entry of
  // the mmm or mm2 list is the same entry of the state list for every pair. The cluster mappings assume it, and so
  // does the cache, whose key is the elements in the order of the state list. Only the state list is built for every
  // pair, the other two orders are taken from the first pair.
  const std::pair<size_t, size_t> first_pair{0, reference_config.GetFirstNeighborsAdjacencyList()[0][0]};
  const std::pair<size_t, size_t> reversed_first_pair{first_pair.second, first_pair.first};
  const auto first_state_list = GetSortedLatticeVectorStateOfPair(reference_config, first_pair);
  const auto positions_in_state_list = [&first_state_list](const std::vector<cfg::Lattice> &lattice_list) {
    std::vector<size_t> positions;
    for (const auto &lattice: lattice_list) {
      const auto it = std::find_if(first_state_list.begin(), first_state_list.end(),
                                   [&lattice](const cfg::Lattice &state_lattice) {
                                     return state_lattice.GetId() == lattice.GetId();
                                   });
      positions.push_back(static_cast<size_t>(it - first_state_list.begin()));
    }
    return positions;
  };
  mmm_positions_ = positions_in_state_list(GetSymmetricallySortedLatticeVectorMMM(reference_config, first_pair));
  mm2_positions_ = positions_in_state_list(GetSymmetricallySortedLatticeVectorMM2(reference_config, first_pair));
  // The frame of a jump is fixed up to a half turn about the jump direction: which of the two neighbours across the
  // jump direction becomes the second axis depends on the order of the neighbour list. The frame of the reversed
  // pair can therefore be turned by half a turn against the frame of the pair, and its mm2 list then names the same
  // sites in another order. The mm2 encoding does not tell the two apart, because it only uses the distance of a site
  // from the centre and its coordinate along the jump. One of the two orders is used for every pair, so that the
  // model is a function of the elements in the order of the state list alone, to the last bit.
  mm2_backward_positions_ =
      positions_in_state_list(GetSymmetricallySortedLatticeVectorMM2(reference_config, reversed_first_pair));
  // Entry k of the state list after half a turn about the jump direction, the x axis of the frame.
  std::vector<size_t> half_turn(first_state_list.size(), first_state_list.size());
  for (size_t k = 0; k < first_state_list.size(); ++k) {
    const auto &position = first_state_list[k].GetRelativePosition();
    for (size_t l = 0; l < first_state_list.size(); ++l) {
      const auto &image = first_state_list[l].GetRelativePosition();
      if (std::abs(image[0] - position[0]) < kEpsilon && std::abs(image[1] - (1 - position[1])) < kEpsilon &&
          std::abs(image[2] - (1 - position[2])) < kEpsilon) {
        half_turn[k] = l;
        break;
      }
    }
  }
  // Entry n of the state list of the reversed pair is entry reversed_positions[n] of the state list of the pair, or
  // the half turn of that entry.
  const auto reversed_positions =
      positions_in_state_list(GetSortedLatticeVectorStateOfPair(reference_config, reversed_first_pair));
  std::vector<size_t> reversed_positions_turned;
  for (const auto position: reversed_positions) {
    reversed_positions_turned.push_back(position < half_turn.size() ? half_turn[position] : half_turn.size());
  }
  const auto is_arranged_like_first_pair = [](const std::vector<size_t> &lattice_id_vector,
                                              const std::vector<size_t> &lattice_id_vector_state,
                                              const std::vector<size_t> &positions) {
    if (lattice_id_vector.size() != positions.size()) {
      return false;
    }
    for (size_t index = 0; index < positions.size(); ++index) {
      if (positions[index] >= lattice_id_vector_state.size() ||
          lattice_id_vector[index] != lattice_id_vector_state[positions[index]]) {
        return false;
      }
    }
    return true;
  };
  bool is_consistent = true;
  // Every site of the three orders has to be a site of the state list: the model reads the state list through them.
  for (const auto *positions: {&mmm_positions_, &mm2_positions_, &mm2_backward_positions_}) {
    for (const auto position: *positions) {
      is_consistent = is_consistent && position < first_state_list.size();
    }
  }
#pragma omp parallel for default(none) shared(reference_config, first_state_list) reduction(&&: is_consistent)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    const auto &neighbor_list = reference_config.GetFirstNeighborsAdjacencyList()[i];
    for (size_t neighbor_slot = 0; neighbor_slot < neighbor_list.size(); ++neighbor_slot) {
      const auto j = neighbor_list[neighbor_slot];
      neighbor_slot_lookup_[i][j] = neighbor_slot;
      const auto flat_idx = i * constants::kNumFirstNearestNeighbors + neighbor_slot;
      const auto sorted_lattice_vector = GetSortedLatticeVectorStateOfPair(reference_config, {i, j});
      // The sites have to lie where those of the first pair lie.
      bool is_like_first_pair = sorted_lattice_vector.size() == first_state_list.size();
      for (size_t k = 0; is_like_first_pair && k < sorted_lattice_vector.size(); ++k) {
        const auto &position = sorted_lattice_vector[k].GetRelativePosition();
        const auto &first_position = first_state_list[k].GetRelativePosition();
        is_like_first_pair = std::abs(position[0] - first_position[0]) < kEpsilon &&
            std::abs(position[1] - first_position[1]) < kEpsilon &&
            std::abs(position[2] - first_position[2]) < kEpsilon;
      }
      is_consistent = is_consistent && is_like_first_pair;
      std::vector<size_t> lattice_id_vector_state;
      lattice_id_vector_state.reserve(sorted_lattice_vector.size());
      std::transform(
          sorted_lattice_vector.begin(),
          sorted_lattice_vector.end(),
          std::back_inserter(lattice_id_vector_state),
          [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_state_flat_[flat_idx] = std::move(lattice_id_vector_state);
    }
  }
  // The frame of the reversed pair has to be one of the two described above.
#pragma omp parallel for default(none) shared(reference_config, is_arranged_like_first_pair, reversed_positions, reversed_positions_turned) reduction(&&: is_consistent)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    const auto &neighbor_list = reference_config.GetFirstNeighborsAdjacencyList()[i];
    for (size_t neighbor_slot = 0; neighbor_slot < neighbor_list.size(); ++neighbor_slot) {
      const auto j = neighbor_list[neighbor_slot];
      const auto reversed_slot = neighbor_slot_lookup_[j].find(i);
      if (reversed_slot == neighbor_slot_lookup_[j].end()) {
        is_consistent = false;
        continue;
      }
      const auto &lattice_id_vector_state =
          site_bond_cluster_state_flat_[i * constants::kNumFirstNearestNeighbors + neighbor_slot];
      const auto &lattice_id_vector_reversed =
          site_bond_cluster_state_flat_[j * constants::kNumFirstNearestNeighbors + reversed_slot->second];
      is_consistent = is_consistent &&
          (is_arranged_like_first_pair(lattice_id_vector_reversed, lattice_id_vector_state, reversed_positions) ||
           is_arranged_like_first_pair(lattice_id_vector_reversed, lattice_id_vector_state,
                                       reversed_positions_turned));
    }
  }
  if (!is_consistent) {
    throw std::runtime_error("The neighbourhoods of the jumps are not arranged alike in this configuration. The "
                             "barrier model needs a complete FCC lattice in a box of at least four unit cells along "
                             "every direction.");
  }
}
VacancyMigrationPredictorQuartic::~VacancyMigrationPredictorQuartic() = default;

std::pair<double, double> VacancyMigrationPredictorQuartic::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
double VacancyMigrationPredictorQuartic::GetDe(const cfg::Config &config,
                                               const std::pair<size_t,
                                                               size_t> &lattice_id_jump_pair) const {
  const auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto flat_idx = GetPairFlatIndex(lattice_id_jump_pair);
  const auto &lattice_id_vector = site_bond_cluster_state_flat_.at(flat_idx);

  // Element codes of the neighbourhood before the jump and after it, when the migrating atom sits on the first site
  // and has left the vacancy on the second one.
  const auto migration_code = state_cluster_indexer_.GetElementCode(migration_element);
  const auto vacancy_code = state_cluster_indexer_.GetElementCode(Element(ElementName::X));
  auto &start_codes = GetThreadLocalStartCodesBuffer();
  auto &end_codes = GetThreadLocalEndCodesBuffer();
  start_codes.resize(lattice_id_vector.size());
  end_codes.resize(lattice_id_vector.size());
  for (size_t i = 0; i < lattice_id_vector.size(); ++i) {
    const size_t lattice_id = lattice_id_vector[i];
    const auto code = state_cluster_indexer_.GetElementCode(config.GetElementAtLatticeId(lattice_id));
    start_codes[i] = code;
    if (lattice_id == lattice_id_jump_pair.first) {
      end_codes[i] = migration_code;
    } else if (lattice_id == lattice_id_jump_pair.second) {
      end_codes[i] = vacancy_code;
    } else {
      end_codes[i] = code;
    }
  }

  // Number of clusters of every kind after the jump minus their number before it. The counts are small integers, so
  // the difference is exact.
  auto &count_changes = GetThreadLocalStartCountsBuffer();
  count_changes.assign(state_cluster_indexer_.Size(), 0);
  AddClusterCountChanges(mapping_state_, state_cluster_indexer_, start_codes, end_codes, count_changes);

  auto &de_encode = GetThreadLocalDeEncodeBuffer();
  de_encode.resize(state_cluster_indexer_.Size());
  const auto &total_bonds = state_cluster_indexer_.GetTotalBonds();
  for (size_t idx = 0; idx < state_cluster_indexer_.Size(); ++idx) {
    de_encode[idx] = static_cast<double>(count_changes[idx]) / total_bonds[idx];
  }

  const Eigen::Map<const Eigen::VectorXd> encode_vec(de_encode.data(), static_cast<Eigen::Index>(de_encode.size()));
  return base_theta_.dot(encode_vec);
}
namespace {
// The dot product of the weights with the one-hot encoding of a neighbourhood. Its sites are those of the state list
// of the pair, taken in the order `positions` gives.
double DotWithEncoding(const OneHotEncoder &encoder,
                       const cfg::Config &config,
                       const std::vector<size_t> &lattice_id_vector_state,
                       const std::vector<size_t> &positions,
                       const std::vector<double> &cluster_weights) {
  auto &element_indices = GetThreadLocalStartCodesBuffer();
  element_indices.resize(positions.size());
  for (size_t i = 0; i < positions.size(); ++i) {
    element_indices[i] = encoder.GetElementIndex(config.GetElementAtLatticeId(lattice_id_vector_state[positions[i]]));
  }
  return encoder.Dot(element_indices.data(), cluster_weights.data());
}
}    // namespace

double VacancyMigrationPredictorQuartic::GetKs(const cfg::Config &config,
                                               const std::pair<size_t,
                                                               size_t> &lattice_id_jump_pair) const {
  const auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);
  const auto &lattice_id_vector_state = site_bond_cluster_state_flat_.at(GetPairFlatIndex(lattice_id_jump_pair));
  // The encoding of the mm2 model is the sum of the encodings seen along the jump and against it.
  const double forward = DotWithEncoding(encoder_mm2_, config, lattice_id_vector_state, mm2_positions_,
                                         element_parameters.cluster_weight_Ks);
  const double backward = DotWithEncoding(encoder_mm2_, config, lattice_id_vector_state, mm2_backward_positions_,
                                          element_parameters.cluster_weight_Ks);
  double logKs = forward + backward - element_parameters.weight_dot_mu_Ks;
  logKs *= element_parameters.sigma_y_Ks;
  logKs += element_parameters.mu_y_Ks;
  return std::exp(logKs);
}
double VacancyMigrationPredictorQuartic::GetD(const cfg::Config &config,
                                              const std::pair<size_t,
                                                              size_t> &lattice_id_jump_pair) const {
  const auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);
  const auto &lattice_id_vector_state = site_bond_cluster_state_flat_.at(GetPairFlatIndex(lattice_id_jump_pair));
  double logD = DotWithEncoding(encoder_mmm_, config, lattice_id_vector_state, mmm_positions_,
                                element_parameters.cluster_weight_D) - element_parameters.weight_dot_mu_D;
  logD *= element_parameters.sigma_y_D;
  logD += element_parameters.mu_y_D;
  return std::exp(logD);
}
double GetQuarticModelBarrier(const double dE, const double D, const double Ks) {
  const auto b = 4 * dE / (D * D * D);
  const auto a = Ks / (4 * D * D);
  const auto c = (9 * b * b - 16 * a * a * D * D) / (32 * a);
  const auto delta = std::sqrt(std::abs(9 * b * b - 32 * a * c));
  const auto Ea = (3 * b + delta) * (3 * b + delta) * (3 * b * b - 16 * a * c + b * delta)
      / std::pow(a, 3) / 2048;
  return Ea;
}
double GetQuarticBarrier(const double dE, const double D, const double Ks) {
  const double model_barrier = GetQuarticModelBarrier(dE, D, Ks);
  const double lowest_barrier = std::max(dE, 0.0);
  // Written as a comparison that is false for a barrier that is not a number, so that such a barrier is returned as
  // it is instead of being replaced by a number.
  return model_barrier <= lowest_barrier ? lowest_barrier : model_barrier;
}
std::pair<double, double> VacancyMigrationPredictorQuartic::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
//     double dE, D, Ks;
// #pragma omp parallel sections default(none) shared(config, lattice_id_jump_pair, dE, D, Ks)
//   {
// #pragma omp section
//     {
//       dE = GetDe(config, lattice_id_jump_pair);
//     }
// #pragma omp section
//     {
//       D = GetD(config, lattice_id_jump_pair);
//     }
// #pragma omp section
//     {
//       Ks = GetKs(config, lattice_id_jump_pair);
//     }
//   }
  const auto dE = GetDe(config, lattice_id_jump_pair);
  const auto D = GetD(config, lattice_id_jump_pair);
  const auto Ks = GetKs(config, lattice_id_jump_pair);
  return {GetQuarticBarrier(dE, D, Ks), dE};
}
double VacancyMigrationPredictorQuartic::GetModelBarrierFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto dE = GetDe(config, lattice_id_jump_pair);
  const auto D = GetD(config, lattice_id_jump_pair);
  const auto Ks = GetKs(config, lattice_id_jump_pair);
  return GetQuarticModelBarrier(dE, D, Ks);
}

size_t VacancyMigrationPredictorQuartic::GetPairFlatIndex(const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto &slot_map = neighbor_slot_lookup_.at(lattice_id_jump_pair.first);
  const auto it = slot_map.find(lattice_id_jump_pair.second);
  if (it == slot_map.end()) {
    throw std::out_of_range("Neighbor not found for lattice pair in GetPairFlatIndex");
  }
  return lattice_id_jump_pair.first * constants::kNumFirstNearestNeighbors + it->second;
}
} // pred
