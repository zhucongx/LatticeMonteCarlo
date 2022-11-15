#ifndef LKMC_LKMC_CFG_INCLUDE_CONFIG_H_
#define LKMC_LKMC_CFG_INCLUDE_CONFIG_H_
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include "Atom.hpp"
#include "Lattice.hpp"
// using Graph = boost::adjacency_list<boost::vecS, boost::vecS>;
namespace cfg {
class Config {
  public:
    /// Constructor
    Config();
    Config(const Matrix_t &basis,
           std::vector<Lattice> lattice_vector,
           std::vector<Atom> atom_vector,
           bool update_neighbor);
    /// Getter
    [[nodiscard]] size_t GetNumAtoms() const;
    [[nodiscard]] const Matrix_t &GetBasis() const;
    [[nodiscard]] const std::vector<Lattice> &GetLatticeVector() const;
    [[nodiscard]] const std::vector<Atom> &GetAtomVector() const;
    [[nodiscard]] const std::vector<std::vector<size_t> > &GetFirstNeighborsAdjacencyList() const;
    [[nodiscard]] const std::vector<std::vector<size_t> > &GetSecondNeighborsAdjacencyList() const;
    [[nodiscard]] const std::vector<std::vector<size_t> > &GetThirdNeighborsAdjacencyList() const;
    [[nodiscard]] std::vector<size_t> GetFirstNeighborsAtomIdVectorOfAtom(size_t atom_id) const;
    [[nodiscard]] std::vector<size_t> GetSecondNeighborsAtomIdVectorOfAtom(size_t atom_id) const;
    [[nodiscard]] size_t GetAtomIdFromLatticeId(size_t lattice_id) const;
    [[nodiscard]] size_t GetLatticeIdFromAtomId(size_t atom_id) const;
    [[nodiscard]] Element GetElementAtAtomId(size_t atom_id) const;
    [[nodiscard]] Element GetElementAtLatticeId(size_t lattice_id) const;
    [[nodiscard]] std::set<Element> GetElementSetWithoutVacancy() const;
    [[nodiscard]] size_t GetStateHash() const;
    /// Modify config
    void AtomJump(const std::pair<size_t, size_t> &atom_id_jump_pair);
    void LatticeJump(const std::pair<size_t, size_t> &lattice_id_jump_pair);
    void ReassignLatticeVector();
    void ChangeAtomElementTypeAtAtom(size_t atom_id, Element element);
    void ChangeAtomElementTypeAtLattice(size_t lattice_id, Element element);
    /// IO
    static Config ReadCfg(const std::string &filename);
    void WriteCfg(const std::string &filename, bool neighbors_info) const;
    static Config ReadMap(const std::string &lattice_filename,
                          const std::string &element_filename,
                          const std::string &map_filename);
    void WriteLattice(const std::string &filename) const;
    void WriteElement(const std::string &filename) const;
    void WriteMap(const std::string &filename) const;
  private:
    /// Private Getter
    [[nodiscard]] const std::unordered_map<size_t, size_t> &GetLatticeToAtomHashmap() const;
    [[nodiscard]] const std::unordered_map<size_t, size_t> &GetAtomToLatticeHashmap() const;
    /// Modify config
    void ConvertRelativeToCartesian();
    void ConvertCartesianToRelative();
    void InitializeNeighborsList(size_t num_atoms);
    void UpdateNeighbors();
    /// Properties
    Matrix_t basis_{};
    std::vector<Lattice> lattice_vector_{};
    std::vector<Atom> atom_vector_{};
    std::unordered_map<size_t, size_t> lattice_to_atom_hashmap_{};
    std::unordered_map<size_t, size_t> atom_to_lattice_hashmap_{};
    // nearest neighbor lists
    std::vector<std::vector<size_t> > first_neighbors_adjacency_list_{};
    std::vector<std::vector<size_t> > second_neighbors_adjacency_list_{};
    std::vector<std::vector<size_t> > third_neighbors_adjacency_list_{};
    // std::vector<std::vector<size_t> > fourth_neighbors_adjacency_list_{};
    // std::vector<std::vector<size_t> > fifth_neighbors_adjacency_list_{};
    // std::vector<std::vector<size_t> > sixth_neighbors_adjacency_list_{};
    // std::vector<std::vector<size_t> > seventh_neighbors_adjacency_list_{};
};
Vector_t GetLatticePairCenter(const Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair);
Matrix_t GetLatticePairRotationMatrix(const Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair);
void RotateLatticeVector(std::vector<Lattice> &lattice_list, const Matrix_t &rotation_matrix);
size_t GetVacancyAtomIndex(const Config &config);
size_t GetVacancyLatticeIndex(const Config &config);
std::unordered_set<size_t> GetNeighborsLatticeIdSetOfJumpPair(
    const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::unordered_set<size_t> GetNeighborsLatticeIdSetOfLatticeId(
    const Config &config, size_t lattice_id);
Config GetNeighborsConfigSetOfJumpPair(
    const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
int FindDistanceLabelBetweenLattice(size_t index1, size_t index2, const Config &config);
Config GenerateFCC(const Factor_t &factors, Element element);
Config GenerateSoluteConfigFromExcitingPure(Config config,
                                            const std::map<Element, size_t> &solute_atom_count);
Config GenerateSoluteConfig(const Factor_t &factors,
                            Element solvent_element,
                            const std::map<Element, size_t> &solute_atom_count);
// Config GenerateClusteredConfigFromExcitingPure(Config config,
//                                             const std::map<Element, size_t> &solute_atom_count);
// Config GenerateClusteredConfig(const Factor_t &factors,
//                             Element solvent_element,
//                             const std::map<Element, size_t> &solute_atom_count);
} // namespace cfg
#endif //LKMC_LKMC_CFG_INCLUDE_CONFIG_H_
