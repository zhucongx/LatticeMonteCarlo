#ifndef LKMC_LKMC_ANSYS_INCLUDE_CLUSTERFINDER_H_
#define LKMC_LKMC_ANSYS_INCLUDE_CLUSTERFINDER_H_
#include <vector>
#include <set>
#include <unordered_set>

#include "Config.h"
#include "EnergyEstimator.h"

namespace ansys {
class ClustersFinder {
  public:
    using ClusterElementNumMap = std::vector<std::pair<std::map<Element, size_t>, double> >;
    ClustersFinder(std::string cfg_filename,
                   Element solvent_atom_type,
                   size_t smallest_cluster_criteria,
                   size_t solvent_bond_criteria,
                   const pred::EnergyEstimator &energy_estimator);

    ClusterElementNumMap FindClustersAndOutput();

  private:
    double GetAbsoluteEnergyOfCluster(const std::vector<size_t> &atom_id_list);
    [[nodiscard]] std::unordered_set<size_t> FindSoluteAtomIndexes() const;
    [[nodiscard]] std::vector<std::vector<size_t> > FindAtomListOfClustersBFSHelper(
        std::unordered_set<size_t> unvisited_atoms_id_set) const;

    // Return a 2D array where values of each row representing the number of atoms of different
    // element in one cluster
    [[nodiscard]] std::vector<std::vector<size_t> > FindAtomListOfClusters() const;

    const std::string cfg_filename_;
    const cfg::Config config_;
    cfg::Config solvent_config_;
    const Element solvent_element_;
    const std::set<Element> element_set_;
    const size_t smallest_cluster_criteria_;
    const size_t solvent_bond_criteria_;

    const pred::EnergyEstimator energy_estimator_;
    // double absolute_energy_solvent_config_{};
};
} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_CLUSTERFINDER_H_
