/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/22/23 12:58 PM                                                          *
 **************************************************************************************************/

/*! \file  Config.h
 *  \brief File for the Config class definition.
 */

#ifndef LMC_CONFIG_INCLUDE_CONFIG_H_
#define LMC_CONFIG_INCLUDE_CONFIG_H_

#include <unordered_map>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <variant>
#include <any>
#include "Eigen/Dense"
#include "Constants.hpp"
#include "Element.hpp"
// #include "Atom.hpp"

/*! \brief Class for defining a configuration of atoms and their positions.
 */
class Config {
 public:
  /*! \brief Default constructor for Config.
   */
  Config();

  /*! \brief Constructor for setting up the configuration of atoms and their positions.
   *  \param basis                    : The basis vectors (3x3) of the configuration.
   *  \param relative_position_matrix : The relative position matrix (3，n) of the configuration.
   *  \param atom_vector              : The atom vector of the configuration，n atoms in total.
   */
  Config(Eigen::Matrix3d basis, Eigen::Matrix3Xd relative_position_matrix, std::vector<Element> atom_vector);

  /*! \brief Default destructor for Config.
   */
  virtual ~Config();

  /*! \brief Query for the number of atoms in the configuration.
   *  \return : The number of atoms in the configuration.
   */
  [[nodiscard]] size_t GetNumAtoms() const;

  /*! \brief Query for the number of lattice sites in the configuration.
   *  \return : The number of lattice sites in the configuration.
   */
  [[nodiscard]] size_t GetNumSites() const;

  /*! \brief Query for the basis vectors of the configuration.
   *  \return : The basis vectors of the configuration.
   */
  [[nodiscard]] const Eigen::Matrix3d &GetBasis() const;

  /*! \brief Query for the lists of neighbor sites of each site in the configuration.
   *  \return : The lists of neighbor sites of each site in the configuration.
   */
  [[nodiscard]] const std::vector<std::vector<std::vector<size_t> > > &GetNeighborLists() const;

  /*! \brief Query for the order of distance between two lattice.
   *  \param lattice_id1 : The lattice id of the first lattice.
   *  \param lattice_id2 : The lattice id of the second lattice.
   *  \return            : The the order of distance between the two lattice. 0 means the same lattice or not neighbors.
   */
  [[nodiscard]] size_t GetDistanceOrder(size_t lattice_id1, size_t lattice_id2) const;

  /*! \brief Query for the cartesian position of an lattice site.
   *  \param lattice_id : The lattice id of the lattice site.
   *  \return           : The cartesian position of the lattice site.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetRelativePositionOfSite(size_t lattice_id) const;

  /*! \brief Query for the cartesian position of an lattice site.
   *  \param lattice_id : The lattice id of the lattice site.
   *  \return           : The cartesian position of the lattice site.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetCartesianPositionOfSite(size_t lattice_id) const;

  /*! \brief Query for the relative position of an atom.
   *  \param atom_id : The atom id of the atom.
   *  \return : The relative position of the atom.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetRelativePositionOfAtom(size_t atom_id) const;

  /*! \brief Query for the cartesian position of an atom.
   *  \param atom_id : The atom id of the atom.
   *  \return        : The cartesian position of the atom.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetCartesianPositionOfAtom(size_t atom_id) const;

  /*! \brief Set the periodic boundary condition of the configuration.
   *  \param periodic_boundary_condition : The periodic boundary condition of the configuration.
   */
  void SetPeriodicBoundaryCondition(const std::array<bool, 3> &periodic_boundary_condition);

  /*! \brief Set the element type at the atom with given atom id.
   *  \param atom_id      : The atom id of the atom.
   *  \param element_type : The new element type.
   */
  void ChangeAtomElementTypeAtAtomId(size_t atom_id, Element element_type);

  /*! \brief Modify the atom configuration.
   *  \param atom_id_jump_pair : The pair of atom ids to modify the configuration.
   */
  void AtomJump(const std::pair<size_t, size_t> &atom_id_jump_pair);

  /*! \brief Update the neighbor list of the configuration with the given cutoffs.
   *  \param cutoffs : The cutoffs to update the neighbor list.
   */
  void UpdateNeighborList(std::vector<double> cutoffs);

  /*! \brief Read the configuration from a lattice file, element file and map file.
   *  \param lattice_filename : The name of the lattice file.
   *  \param element_filename : The name of the element file.
   *  \param map_filename     : The name of the map file.
   *  \return                 : The configuration read from the file.
   */
  static Config ReadMap(const std::string &lattice_filename,
                        const std::string &element_filename,
                        const std::string &map_filename);

  /*! \brief Read the configuration from a CFG file.
   *  \param filename : The name of the CFG file.
   *  \return         : The configuration read from the file.
   */
  static Config ReadCfg(const std::string &filename);

  /*! \brief Read the configuration from a POSCAR file.
   *  \param filename : The name of the POSCAR file.
   *  \return         : The configuration read from the file.
   */
  static Config ReadPoscar(const std::string &filename);

  /*! \brief Write the configuration to a file.
   *  \param filename : The name of the file to write the configuration to.
   */
  void WriteConfig(const std::string &filename) const;

  /*! \brief Write the extended configuration to a file. Check the extended CFG format
   *         at <http://li.mit.edu/Archive/Graphics/A/#extended_CFG>.
   *  \param filename        : The name of the file to write the extended configuration to.
   *  \param auxiliary_lists : The auxiliary lists to write to the file.
   */
  void WriteConfigExtended(
      const std::string &filename,
      const std::map<std::string, std::vector<double>> &auxiliary_lists) const;

  using VectorVariant = std::variant<std::vector<int>, std::vector<size_t>, std::vector<double>,
                                     std::vector<std::string>, std::vector<Eigen::Vector3d> >;
  using ValueVariant = std::variant<int, double, std::string>;

  /*! \brief Write the extended XXY to a file. Check the extended XYZ format
   *         at <https://web.archive.org/web/20190811094343/https://libatoms.github.io/QUIP/io.html#extendedxyz>
   *         and <https://www.ovito.org/docs/current/reference/file_formats/input/xyz.html#file-formats-input-xyz>
   *  \param filename        : The name of the file to write the extended configuration to.
   *  \param auxiliary_lists : The lists of atom information .
   *  \param global_list     : The list of global information.
   */
  void WriteXyzExtended(
      const std::string &filename,
      const std::map<std::string, VectorVariant> &auxiliary_lists,
      const std::map<std::string, ValueVariant> &global_list) const;

 private:
  /*! \brief Sort lattice sites by the positions (x, y, z)
   */
  void ReassignLattice();

  /*! \brief Query for the relative distance vector between two lattice.
   *  \param lattice_id1 : The lattice id of the first lattice.
   *  \param lattice_id2 : The lattice id of the second lattice.
   *  \return            : The relative distance vector between the two lattice.
   */
  [[nodiscard]] Eigen::Vector3d GetRelativeDistanceVectorLattice(size_t lattice_id1, size_t lattice_id2) const;

  /// The periodic boundary condition status in all three directions.
  std::array<bool, 3> periodic_boundary_condition_{true, true, true};

  /// The basis matrix of the configuration.
  Eigen::Matrix3d basis_{};

  /// The relative position matrix of the configuration.
  Eigen::Matrix3Xd relative_position_matrix_{};

  /// The cartesian position matrix of the configuration.
  Eigen::Matrix3Xd cartesian_position_matrix_{};

  /// The vector of atoms in the configuration.
  std::vector<Element> atom_vector_{};

  /// Mapping from lattice points to atom ids.
  std::unordered_map<size_t, size_t> lattice_to_atom_hashmap_{};

  /// Mapping from atom ids to lattice points.
  std::unordered_map<size_t, size_t> atom_to_lattice_hashmap_{};

  /// Nearest neighbor lists. The first key is the cutoff distance between the two lattice sites in Angstrom.
  std::vector<std::vector<std::vector<size_t> > > neighbor_lists_{};

  /// Cutoffs for the neighbor lists.
  std::vector<double> cutoffs_{};

  /// Number of cells in each direction.
  Eigen::Vector3i num_cells_{};

  /// The cells of the configuration. Each cell is a vector of lattice ids.
  std::vector<std::vector<size_t>> cells_{};
};

#endif //LMC_CONFIG_INCLUDE_CONFIG_H_
