#ifndef LMC_LMC_PRED_TEST_SYNTHETICMODEL_HPP_
#define LMC_LMC_PRED_TEST_SYNTHETICMODEL_HPP_
// Helpers shared by the tests that need a configuration and a coefficient file.
#include "Config.h"
#include "EnergyUtility.h"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace pred {
inline const std::set<Element> &TestElementSet() {
  static const std::set<Element> element_set{Element("Al"), Element("Mg"), Element("Zn")};
  return element_set;
}

// 256 sites of Al of which about 30 % are replaced by Mg and Zn, with one vacancy at site 0.
inline cfg::Config MakeTestAlloy() {
  auto config = cfg::GenerateFCC({4, 4, 4}, Element("Al"));
  std::mt19937_64 generator(1);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  for (size_t lattice_id = 1; lattice_id < config.GetNumAtoms(); ++lattice_id) {
    const double draw = uniform(generator);
    if (draw < 0.15) {
      config.SetAtomElementTypeAtLattice(lattice_id, Element("Mg"));
    } else if (draw < 0.30) {
      config.SetAtomElementTypeAtLattice(lattice_id, Element("Zn"));
    }
  }
  config.SetAtomElementTypeAtLattice(0, Element("X"));
  return config;
}

// Returns the layout of a coefficient file for the element set, filled with random numbers.
//
// The numbers have no physical meaning. They are scaled so that the energy change of a jump is of the order of
// 0.01 to 0.1 eV and barriers are about 0.3 eV. Every relation between the predictors that the tests check has to
// hold for any coefficients.
inline nlohmann::json MakeSyntheticModel(const cfg::Config &reference_config, const std::set<Element> &element_set) {
  std::mt19937_64 generator(2024);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  const auto random_vector = [&](const size_t size, const double scale, const double offset) {
    std::vector<double> values(size);
    for (auto &value: values) {
      value = offset + scale * uniform(generator);
    }
    return values;
  };
  const auto random_matrix = [&](const size_t num_rows, const size_t num_columns, const double scale) {
    std::vector<std::vector<double>> rows;
    for (size_t row = 0; row < num_rows; ++row) {
      rows.push_back(random_vector(num_columns, scale, 0.0));
    }
    return rows;
  };

  auto element_set_with_vacancy(element_set);
  element_set_with_vacancy.emplace(ElementName::X);
  const size_t encode_size_mmm =
      GetOneHotEncodeSize(element_set.size(), GetAverageClusterParametersMappingMMM(reference_config));
  const size_t encode_size_mm2 =
      GetOneHotEncodeSize(element_set.size(), GetAverageClusterParametersMappingMM2(reference_config));
  const size_t num_components = 4;

  nlohmann::json model;
  model["Base"]["theta"] = random_vector(InitializeClusterHashMap(element_set_with_vacancy).size(), 0.5, 0.0);
  for (const auto &element: element_set) {
    auto &parameters = model[element.GetString()];
    parameters["mu_x_mmm"] = random_vector(encode_size_mmm, 0.1, 0.0);
    parameters["mu_x_mm2"] = random_vector(encode_size_mm2, 0.1, 0.0);
    parameters["sigma_x_mmm"] = random_vector(encode_size_mmm, 0.2, 1.0);
    parameters["sigma_x_mm2"] = random_vector(encode_size_mm2, 0.2, 1.0);
    parameters["U_mmm"] = random_matrix(num_components, encode_size_mmm, 0.05);
    parameters["U_mm2"] = random_matrix(num_components, encode_size_mm2, 0.05);
    parameters["theta_D"] = random_vector(num_components, 0.5, 0.0);
    parameters["theta_Ks"] = random_vector(num_components, 0.5, 0.0);
    parameters["mu_D"] = 1.0;
    parameters["mu_Ks"] = 1.0;
    parameters["sigma_D"] = 0.05;
    parameters["sigma_Ks"] = 0.05;
  }
  return model;
}

inline std::string WriteModel(const std::filesystem::path &filename, const nlohmann::json &model) {
  std::ofstream ofs(filename);
  ofs << model.dump();
  return filename.string();
}
}    // namespace pred
#endif    //LMC_LMC_PRED_TEST_SYNTHETICMODEL_HPP_
