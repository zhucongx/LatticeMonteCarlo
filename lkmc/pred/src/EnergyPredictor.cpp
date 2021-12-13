#include "EnergyPredictor.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &type_set) {
  size_t type_size = type_set.size();
  std::unordered_map<std::string, std::vector<double> > encode_dict;

  size_t ct1 = 0;
  for (const auto &element: type_set) {
    std::vector<double> element_encode(type_size, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element.GetString()] = element_encode;
    ++ct1;
  }

  size_t num_pairs = type_size * type_size;
  size_t ct2 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      std::vector<double> element_encode(num_pairs, 0);
      element_encode[ct2] = 1.0;
      encode_dict[element1.GetString() + element2.GetString()] = element_encode;
      ++ct2;
    }
  }

  size_t num_triplets = type_size * type_size * type_size;
  size_t ct3 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      for (const auto &element3: type_set) {
        std::vector<double> element_encode(num_triplets, 0);
        element_encode[ct3] = 1.0;
        encode_dict[element1.GetString() + element2.GetString() + element3.GetString()] =
            element_encode;
        ++ct3;
      }
    }
  }
  return encode_dict;
}
std::vector<double> GetOneHotParametersFromMap(
    const std::vector<Element> &encode,
    const std::unordered_map<std::string, std::vector<double> > &one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {

  std::vector<double> res_encode;
  res_encode.reserve(882);

  for (const auto &cluster_vector: cluster_mapping) {
    std::vector<double> sum_of_list(
        static_cast<size_t>(std::pow(num_of_elements, cluster_vector[0].size())), 0);
    for (const auto &cluster: cluster_vector) {
      std::string cluster_type;
      for (auto index: cluster) {
        cluster_type += encode[index].GetString();
      }
      const auto &cluster_one_hot_encode = one_hot_encode_hashmap.at(cluster_type);
      std::transform(sum_of_list.begin(), sum_of_list.end(),
                     cluster_one_hot_encode.begin(),
                     sum_of_list.begin(),
                     std::plus<>());
    }
    auto cluster_vector_size = static_cast<double>( cluster_vector.size());
    std::for_each(sum_of_list.begin(),
                  sum_of_list.end(),
                  [cluster_vector_size](auto &n) { n /= cluster_vector_size; });

    std::move(sum_of_list.begin(), sum_of_list.end(), std::back_inserter(res_encode));
  }
  return res_encode;
}

EnergyPredictor::EnergyPredictor(const std::set<Element> &type_set)
    : type_set_(type_set),
      one_hot_encode_hash_map_(GetOneHotEncodeHashmap(type_set)) {}
EnergyPredictor::~EnergyPredictor() = default;
std::pair<double, double> EnergyPredictor::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {

  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}

} // namespace pred
