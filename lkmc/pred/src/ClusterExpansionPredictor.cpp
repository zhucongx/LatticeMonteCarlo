#include "ClusterExpansionPredictor.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
ClusterExpansionPredictor::ClusterExpansionPredictor(const std::string &predictor_filename,
                                                     const cfg::Config &reference_config,
                                                     const std::set<std::string> &type_set) {
  // std::ifstream ifs(predictor_filename, std::ifstream::in);
  // json all_parameters;
  // ifs >> all_parameters;
  // for (const auto &[element, parameters]: all_parameters.items()) {
  //   if (element == "Bond") {
  //     theta_ = std::vector<double>(parameters.at("theta"));
  //     continue;
  //   }
  //   element_parameters_hashmap_[element] = Element_Parameters{
  //       parameters.at("mu_x"),
  //       parameters.at("transform_matrix"),
  //       parameters.at("theta"),
  //       parameters.at("mean_y")};
  // }
}
ClusterExpansionPredictor::~ClusterExpansionPredictor() = default;

} // namespace pred