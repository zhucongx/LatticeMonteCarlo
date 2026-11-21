#include "VacancyMigrationPredictorQuartic.h"
#include <Eigen/Dense>
#include <utility>
#include <boost/range/combine.hpp>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "Constants.hpp"
#include <stdexcept>
#include <algorithm>

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
    const auto &parameters = it->second;
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
  }
  site_bond_cluster_state_flat_.resize(num_pair_slots_);
  site_bond_cluster_mmm_flat_.resize(num_pair_slots_);
  site_bond_cluster_mm2_flat_.resize(num_pair_slots_);
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    const auto &neighbor_list = reference_config.GetFirstNeighborsAdjacencyList()[i];
    for (size_t neighbor_slot = 0; neighbor_slot < neighbor_list.size(); ++neighbor_slot) {
      const auto j = neighbor_list[neighbor_slot];
      neighbor_slot_lookup_[i][j] = neighbor_slot;
      const auto flat_idx = i * constants::kNumFirstNearestNeighbors + neighbor_slot;
      auto sorted_lattice_vector =
          GetSortedLatticeVectorStateOfPair(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_state;
      std::transform(
          sorted_lattice_vector.begin(),
          sorted_lattice_vector.end(),
          std::back_inserter(lattice_id_vector_state),
          [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_state_flat_[flat_idx] = lattice_id_vector_state;
      auto sorted_lattice_vector_mmm = GetSymmetricallySortedLatticeVectorMMM(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mmm;
      std::transform(
          sorted_lattice_vector_mmm.begin(),
          sorted_lattice_vector_mmm.end(),
          std::back_inserter(lattice_id_vector_mmm),
          [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mmm_flat_[flat_idx] = lattice_id_vector_mmm;
      auto sorted_lattice_vector_mm2 = GetSymmetricallySortedLatticeVectorMM2(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mm2;
      std::transform(
          sorted_lattice_vector_mm2.begin(),
          sorted_lattice_vector_mm2.end(),
          std::back_inserter(lattice_id_vector_mm2),
          [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mm2_flat_[flat_idx] = lattice_id_vector_mm2;
    }
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
// Writes the one-hot encoding of the elements on the given sites to `out_encode`.
void EncodeNeighbourhood(const OneHotEncoder &encoder,
                         const cfg::Config &config,
                         const std::vector<size_t> &lattice_id_vector,
                         std::vector<double> &out_encode) {
  auto &element_indices = GetThreadLocalStartCodesBuffer();
  element_indices.resize(lattice_id_vector.size());
  for (size_t i = 0; i < lattice_id_vector.size(); ++i) {
    element_indices[i] = encoder.GetElementIndex(config.GetElementAtLatticeId(lattice_id_vector[i]));
  }
  out_encode.resize(encoder.Size());
  encoder.Encode(element_indices.data(), out_encode.data());
}
}    // namespace

double VacancyMigrationPredictorQuartic::GetKs(const cfg::Config &config,
                                               const std::pair<size_t,
                                                               size_t> &lattice_id_jump_pair) const {
  const auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);

  auto &encode_mm2_forward = GetThreadLocalEncodeMM2ForwardBuffer();
  EncodeNeighbourhood(encoder_mm2_, config, site_bond_cluster_mm2_flat_.at(GetPairFlatIndex(lattice_id_jump_pair)),
                      encode_mm2_forward);
  auto &encode_mm2_backward = GetThreadLocalEncodeMM2BackwardBuffer();
  EncodeNeighbourhood(
      encoder_mm2_, config,
      site_bond_cluster_mm2_flat_.at(GetPairFlatIndex({lattice_id_jump_pair.second, lattice_id_jump_pair.first})),
      encode_mm2_backward);

  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);
  // Vectorized normalization using Eigen.
  Eigen::Map<Eigen::VectorXd> encode_vec(encode_mm2_forward.data(), static_cast<Eigen::Index>(encode_mm2_forward.size()));
  const Eigen::Map<const Eigen::VectorXd> backward_vec(encode_mm2_backward.data(),  static_cast<Eigen::Index>(encode_mm2_backward.size()));
  encode_vec += backward_vec;
  encode_vec -= element_parameters.mu_x_mm2;
  encode_vec = encode_vec.cwiseQuotient(element_parameters.sigma_x_mm2);

  const auto &U_mat = element_parameters.U_mm2;
  double logKs = (U_mat * encode_vec).dot(element_parameters.theta_Ks);
  logKs *= element_parameters.sigma_y_Ks;
  logKs += element_parameters.mu_y_Ks;
  return std::exp(logKs);
}
double VacancyMigrationPredictorQuartic::GetD(const cfg::Config &config,
                                              const std::pair<size_t,
                                                              size_t> &lattice_id_jump_pair) const {
  const auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto &encode_mmm = GetThreadLocalEncodeMMMBuffer();
  EncodeNeighbourhood(encoder_mmm_, config, site_bond_cluster_mmm_flat_.at(GetPairFlatIndex(lattice_id_jump_pair)),
                      encode_mmm);

  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  // Vectorized normalization using Eigen.
  Eigen::Map<Eigen::VectorXd> encode_vec(encode_mmm.data(), static_cast<Eigen::Index>(encode_mmm.size()));
  encode_vec -= element_parameters.mu_x_mmm;
  encode_vec = encode_vec.cwiseQuotient(element_parameters.sigma_x_mmm);

  const auto &U_mat = element_parameters.U_mmm;
  double logD = (U_mat * encode_vec).dot(element_parameters.theta_D);
  logD *= element_parameters.sigma_y_D;
  logD += element_parameters.mu_y_D;
  return std::exp(logD);
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
  const auto b = 4 * dE / (D * D * D);
  const auto a = Ks / (4 * D * D);
  const auto c = (9 * b * b - 16 * a * a * D * D) / (32 * a);
  const auto delta = std::sqrt(std::abs(9 * b * b - 32 * a * c));
  const auto Ea = (3 * b + delta) * (3 * b + delta) * (3 * b * b - 16 * a * c + b * delta)
      / std::pow(a, 3) / 2048;
  return {Ea, dE};
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
