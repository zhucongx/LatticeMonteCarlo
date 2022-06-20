#ifndef LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
#define LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
#include "EnergyEstimator.h"
#include "ClustersFinder.h"
namespace ansys {

class Iterator {
  public:
    Iterator(unsigned long long int initial_number,
             unsigned long long int increment_number,
             Element solvent_element,
             std::set<Element> type_set,
             size_t smallest_cluster_criteria,
             size_t solvent_bond_criteria,
             const std::string &predictor_filename);
    virtual ~Iterator();
    void SerialRunCluster() const;
    void SerialRunReformat() const;
  private:
    const unsigned long long initial_number_;
    const unsigned long long increment_number_;
    unsigned long long final_number_;
    Element solvent_element_;
    size_t smallest_cluster_criteria_;
    size_t solvent_bond_criteria_;
    std::unordered_map<unsigned long long, double> filename_time_hashset_{};
    const pred::EnergyEstimator energy_estimator_;
};

} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
