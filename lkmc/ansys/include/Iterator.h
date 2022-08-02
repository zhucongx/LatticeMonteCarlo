#ifndef LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
#define LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
#include "TotalEnergyPredictor.h"
#include "ClustersFinder.h"
namespace ansys {

class Iterator {
  public:
    Iterator(unsigned long long int initial_steps,
             unsigned long long int increment_steps,
             Element solvent_element,
             std::set<Element> element_set,
             size_t smallest_cluster_criteria,
             // size_t solvent_bond_criteria,
             const std::string &predictor_filename);
    virtual ~Iterator();
    void RunCluster() const;
    void RunReformat() const;
  private:
    const unsigned long long initial_steps_;
    const unsigned long long increment_steps_;
    unsigned long long final_number_;
    Element solvent_element_;
    size_t smallest_cluster_criteria_;
    // size_t solvent_bond_criteria_;
    std::unordered_map<unsigned long long, double> filename_time_hashset_{};
    const pred::TotalEnergyPredictor energy_estimator_;
};

} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
