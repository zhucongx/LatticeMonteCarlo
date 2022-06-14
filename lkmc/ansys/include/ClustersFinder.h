#ifndef LKMC_LKMC_ANSYS_INCLUDE_CLUSTERFINDER_H_
#define LKMC_LKMC_ANSYS_INCLUDE_CLUSTERFINDER_H_
#include <vector>
#include <set>
#include <unordered_set>

#include "Config.h"
namespace ansys {
class ClustersFinder {
  public:
    using ClusterElementNumMap = std::vector<std::map<Element, size_t> >;
    ClustersFinder(std::string cfg_filename,
                   Element solvent_atom_type,
                   size_t smallest_cluster_criteria,
                   size_t solvent_bond_criteria);

    ClusterElementNumMap FindClustersAndOutput();

    // static void PrintLog(const unsigned long long int &file_index,
    //                      double time,
    //                      const ClusterElementNumMap &found_data);
  private:
    [[nodiscard]] std::unordered_set<size_t> FindSoluteAtomIndexes() const;
    [[nodiscard]] std::vector<std::vector<size_t> > FindAtomListOfClustersBFSHelper(
        std::unordered_set<size_t> unvisited_atoms_id_set) const;

    // Return a 2D array where values of each row representing the number of atoms of different
    // element in one cluster
    [[nodiscard]] std::vector<std::vector<size_t> > FindAtomListOfClusters() const;

    const std::string cfg_filename_;
    const cfg::Config config_;
    const Element solvent_element_;
    const std::set<Element> element_set_;
    const size_t smallest_cluster_criteria_{};
    const size_t solvent_bond_criteria_{};
};
} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_CLUSTERFINDER_H_
