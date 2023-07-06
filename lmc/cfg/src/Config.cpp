/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/6/23 3:01 PM                                                            *
 **************************************************************************************************/

/*! \file  Config.cpp
 *  \brief File for the Config class implementation.
 */

#include "Config.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

Config::Config() = default;

Config::Config(Eigen::Matrix3d basis,
               Eigen::Matrix3Xd relative_position_matrix,
               std::vector<Element> atom_vector)
    : basis_(std::move(basis)),
      relative_position_matrix_(std::move(relative_position_matrix)),
      atom_vector_(std::move(atom_vector)) {
  // Check that the number of atoms is the same in the basis vectors and the relative position matrix
  if (atom_vector_.size() != static_cast<size_t>(relative_position_matrix_.cols())) {
    throw std::runtime_error("Lattice vector and atom vector size do not match, lattice_vector.size() = " +
        std::to_string(relative_position_matrix_.cols()) + ", atom_vector_.size() = " +
        std::to_string(atom_vector_.size()));
  }
  cartesian_position_matrix_ = basis_ * relative_position_matrix_;
  for (size_t id = 0; id < atom_vector_.size(); ++id) {
    // id here is also lattice_id
    lattice_to_atom_hashmap_.emplace(id, id);
    atom_to_lattice_hashmap_.emplace(id, id);
  }
}

Config::~Config() = default;

size_t Config::GetNumAtoms() const {
  return atom_vector_.size();
}

size_t Config::GetNumSites() const {
  return static_cast<size_t>(relative_position_matrix_.cols());
}

const Eigen::Matrix3d &Config::GetBasis() const {
  return basis_;
}

const std::vector<std::vector<std::vector<size_t> > > &Config::GetNeighborLists() const {
  return neighbor_lists_;
}

Eigen::Vector3d Config::GetRelativePositionOfAtom(size_t atom_id) const {
  return relative_position_matrix_.col(static_cast<int>(atom_to_lattice_hashmap_.at(atom_id)));
}

Eigen::Vector3d Config::GetCartesianPositionOfAtom(size_t atom_id) const {
  return cartesian_position_matrix_.col(static_cast<int>(atom_to_lattice_hashmap_.at(atom_id)));
}

void Config::SetPeriodicBoundaryCondition(const std::array<bool, 3> &periodic_boundary_condition) {
  periodic_boundary_condition_ = periodic_boundary_condition;
}

void Config::ChangeAtomElementTypeAtAtomId(size_t atom_id, Element element_type) {
  atom_vector_.at(atom_id) = Element(element_type);
}

void Config::AtomJump(const std::pair<size_t, size_t> &atom_id_jump_pair) {
  const auto [atom_id_lhs, atom_id_rhs] = atom_id_jump_pair;
  const auto lattice_id_lhs = atom_to_lattice_hashmap_.at(atom_id_lhs);
  const auto lattice_id_rhs = atom_to_lattice_hashmap_.at(atom_id_rhs);

  atom_to_lattice_hashmap_.at(atom_id_lhs) = lattice_id_rhs;
  atom_to_lattice_hashmap_.at(atom_id_rhs) = lattice_id_lhs;
  lattice_to_atom_hashmap_.at(lattice_id_lhs) = atom_id_rhs;
  lattice_to_atom_hashmap_.at(lattice_id_rhs) = atom_id_lhs;
}

void Config::ReassignLattice() {

  std::vector<std::pair<Eigen::Vector3d, size_t>> new_lattice_id_vector(GetNumSites());
  for (size_t i = 0; i < GetNumSites(); ++i) {
    new_lattice_id_vector[i] = std::make_pair(relative_position_matrix_.col(static_cast<int>(i)), i);
  }
  std::sort(new_lattice_id_vector.begin(), new_lattice_id_vector.end(),
            [](const auto &lhs, const auto &rhs) -> bool {
              const double x_diff = lhs.first[0] - rhs.first[0];
              if (x_diff < -constants::kEpsilon) { return true; }
              if (x_diff > constants::kEpsilon) { return false; }
              const double y_diff = lhs.first[1] - rhs.first[1];
              if (y_diff < -constants::kEpsilon) { return true; }
              if (y_diff > constants::kEpsilon) { return false; }
              return lhs.first[2] < rhs.first[2] - constants::kEpsilon;
            });

  std::unordered_map<size_t, size_t> old_lattice_id_to_new;
  for (size_t lattice_id = 0; lattice_id < GetNumAtoms(); ++lattice_id) {
    old_lattice_id_to_new.emplace(new_lattice_id_vector.at(lattice_id).second, lattice_id);
  }
  std::unordered_map<size_t, size_t> new_lattice_to_atom_hashmap, new_atom_to_lattice_hashmap;
  for (size_t atom_id = 0; atom_id < GetNumAtoms(); ++atom_id) {
    auto old_lattice_id = atom_to_lattice_hashmap_.at(atom_id);
    auto new_lattice_id = old_lattice_id_to_new.at(old_lattice_id);
    new_lattice_to_atom_hashmap.emplace(new_lattice_id, atom_id);
    new_atom_to_lattice_hashmap.emplace(atom_id, new_lattice_id);
  }
  // Copy sorted columns back into matrix
  for (size_t i = 0; i < GetNumSites(); ++i) {
    relative_position_matrix_.col(static_cast<int>(i)) = new_lattice_id_vector[i].first;
  }
  cartesian_position_matrix_ = basis_ * relative_position_matrix_;
  lattice_to_atom_hashmap_ = new_lattice_to_atom_hashmap;
  atom_to_lattice_hashmap_ = new_atom_to_lattice_hashmap;
}

Eigen::Vector3d Config::GetRelativeDistanceVectorLattice(size_t lattice_id1, size_t lattice_id2) const {
  Eigen::Vector3d relative_distance_vector = relative_position_matrix_.col(static_cast<int>(lattice_id2))
      - relative_position_matrix_.col(static_cast<int>(lattice_id1));
  // periodic boundary conditions
  for (const size_t kDim : std::vector<size_t>{0, 1, 2}) {
    if (periodic_boundary_condition_[kDim]) {
      while (relative_distance_vector[static_cast<int>(kDim)] >= 0.5) {
        relative_distance_vector[static_cast<int>(kDim)] -= 1;
      }
      while (relative_distance_vector[static_cast<int>(kDim)] < -0.5) {
        relative_distance_vector[static_cast<int>(kDim)] += 1;
      }
    }
  }
  return relative_distance_vector;
}
size_t Config::GetDistanceOrder(size_t lattice_id1, size_t lattice_id2) const {
  if (lattice_id1 == lattice_id2) { return 0; } // same lattice
  Eigen::Vector3d relative_distance_vector = GetRelativeDistanceVectorLattice(lattice_id1, lattice_id2);
  double cartesian_distance = (basis_ * relative_distance_vector).norm();
  auto upper = std::upper_bound(cutoffs_.begin(), cutoffs_.end(), cartesian_distance);
  if (upper == cutoffs_.end()) {
    return std::numeric_limits<double>::infinity(); // distance is larger than the largest cutoff
  }
  return static_cast<size_t>(1 + std::distance(cutoffs_.begin(), upper));
}
void Config::UpdateNeighborList(std::vector<double> cutoffs) {
  cutoffs_ = std::move(cutoffs);
  std::sort(cutoffs_.begin(), cutoffs_.end());
  std::vector<double> cutoffs_squared(cutoffs_.size());
  std::transform(cutoffs_.begin(),
                 cutoffs_.end(),
                 cutoffs_squared.begin(),
                 [](double cutoff) { return cutoff * cutoff; });

  std::vector<std::vector<size_t>> neighbors_list(GetNumSites());
  neighbor_lists_ = {cutoffs_.size(), neighbors_list};

  const auto basis_inverse = basis_.inverse().transpose();
  double cutoff_pbc = std::numeric_limits<double>::infinity();
  for (const auto kDim : {0, 1, 2}) {
    if (periodic_boundary_condition_[static_cast<size_t>(kDim)]) {
      cutoff_pbc = std::min(cutoff_pbc, 0.5 / basis_inverse.col(kDim).norm());
    }
  }
  if (cutoff_pbc <= cutoffs_.back()) {
    throw std::runtime_error(
        "The cutoff is larger than the maximum cutoff allowed due to periodic boundary conditions, "
        "cutoff_pbc = " + std::to_string(cutoff_pbc) + ", cutoff_input = " + std::to_string(cutoffs_.back()));
  }
  // Calculate the number of cells in each dimension, at least 3
  num_cells_ = ((basis_.colwise().norm().array()) / (0. + cutoffs_.back())).floor().cast<int>();
  num_cells_ = num_cells_.cwiseMax(Eigen::Vector3i::Constant(3));

  // Create cells
  cells_ = std::vector<std::vector<size_t>>(static_cast<size_t>(num_cells_.prod()));
  for (size_t lattice_id = 0; lattice_id < GetNumSites(); ++lattice_id) {
    const Eigen::Vector3d relative_position = relative_position_matrix_.col(static_cast<int>(lattice_id));
    Eigen::Vector3i cell_pos = (num_cells_.cast<double>().array() * relative_position.array()).floor().cast<int>();
    int cell_idx = (cell_pos(0) * num_cells_(1) + cell_pos(1)) * num_cells_(2) + cell_pos(2);
    cells_.at(static_cast<size_t>(cell_idx)).push_back(lattice_id);
  }

  // // Check if the max distance between two atoms in the neighboring cells is greater than the cutoff
  // Eigen::Matrix3d cell_basis = basis_.array().colwise() / num_cells_.cast<double>().array();
  // double cutoff_cell = 1 / (cell_basis.inverse().colwise().norm().maxCoeff());
  // if (cutoff_cell <= cutoffs_.back()) {
  //   throw std::runtime_error(
  //       "The cutoff is larger than the maximum cutoff allowed due to non-orthogonal cells "
  //       "cutoff_cell = " + std::to_string(cutoff_cell) + ", cutoff_input = " + std::to_string(cutoffs_.back()));
  // }

  // Create neighbor list, iterate over each cell and find neighboring points
  for (int i = 0; i < num_cells_(0); ++i) {
    for (int j = 0; j < num_cells_(1); ++j) {
      for (int k = 0; k < num_cells_(2); ++k) {
        int cell_idx = (i * num_cells_(1) + j) * num_cells_(2) + k;
        auto &cell = cells_[static_cast<size_t>(cell_idx)];
        // Check neighboring cells, taking into account periodic boundaries
        for (int di = -1; di <= 1; ++di) {
          for (int dj = -1; dj <= 1; ++dj) {
            for (int dk = -1; dk <= 1; ++dk) {
              int ni = (i + di + num_cells_(0)) % num_cells_(0);
              int nj = (j + dj + num_cells_(1)) % num_cells_(1);
              int nk = (k + dk + num_cells_(2)) % num_cells_(2);
              int neighbor_cell_idx = (ni * num_cells_(1) + nj) * num_cells_(2) + nk;
              auto &neighbor_cell = cells_.at(static_cast<size_t>(neighbor_cell_idx));
              // For each point in the cell, check if it's close to any point in the neighboring cell
              for (size_t lattice_id1 : cell) {
                for (size_t lattice_id2 : neighbor_cell) {
                  // Make sure we're not comparing a point to itself, and don't double-count pairs within the same cell
                  if (lattice_id2 >= lattice_id1) { continue; }
                  // Calculate distance
                  Eigen::Vector3d relative_distance_vector = GetRelativeDistanceVectorLattice(lattice_id1, lattice_id2);
                  double cartesian_distance_squared = (basis_ * relative_distance_vector).squaredNorm();
                  // If the distance is less than the cutoff, the points are bonded
                  for (size_t cutoff_squared_id = 0; cutoff_squared_id < cutoffs_squared.size(); ++cutoff_squared_id) {
                    if (cartesian_distance_squared < cutoffs_squared.at(cutoff_squared_id)) {
                      neighbor_lists_.at(cutoff_squared_id).at(lattice_id1).push_back(lattice_id2);
                      neighbor_lists_.at(cutoff_squared_id).at(lattice_id2).push_back(lattice_id1);
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
Config Config::ReadMap(const std::string &lattice_filename,
                       const std::string &element_filename,
                       const std::string &map_filename) {
  std::ifstream ifs_lattice(lattice_filename, std::ifstream::in);
  if (!ifs_lattice) {
    throw std::runtime_error("Cannot open " + lattice_filename);
  }
  size_t num_atoms;
  ifs_lattice >> num_atoms;
  ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  Eigen::Matrix3d basis;
  ifs_lattice >> basis(0, 0) >> basis(0, 1) >> basis(0, 2); // lattice vector a
  ifs_lattice >> basis(1, 0) >> basis(1, 1) >> basis(1, 2); // lattice vector b
  ifs_lattice >> basis(2, 0) >> basis(2, 1) >> basis(2, 2); // lattice vector c
  ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line

  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);
  double position_X, position_Y, position_Z;
  for (size_t lattice_id = 0; lattice_id < num_atoms; ++lattice_id) {
    ifs_lattice >> position_X >> position_Y >> position_Z;
    ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    relative_position_matrix.col(static_cast<int>(lattice_id)) = Eigen::Vector3d{position_X, position_Y, position_Z};
  }

  std::ifstream ifs_element(element_filename, std::ifstream::in);
  if (!ifs_element) {
    throw std::runtime_error("Cannot open " + element_filename);
  }

  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  std::string type;
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id) {
    ifs_element >> type;
    atom_vector.emplace_back(type);
    ifs_element.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  Config config(basis, relative_position_matrix, atom_vector);
  std::ifstream ifs_map(map_filename, std::ifstream::in);
  if (!ifs_map) {
    throw std::runtime_error("Cannot open " + map_filename);
  }
  size_t lattice_id;
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id) {
    ifs_map >> lattice_id;
    config.lattice_to_atom_hashmap_.emplace(lattice_id, atom_id);
    config.atom_to_lattice_hashmap_.emplace(atom_id, lattice_id);
    ifs_map.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  config.ReassignLattice();
  return config;
}
Config Config::ReadCfg(const std::string &filename) {
  std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
  if (!ifs) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  boost::iostreams::filtering_istream fis;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fis.push(boost::iostreams::gzip_decompressor());
  } else if (boost::filesystem::path(filename).extension() == ".bz2") {
    fis.push(boost::iostreams::bzip2_decompressor());
  }
  fis.push(ifs);

  // "Number of particles = %i"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  size_t num_atoms;
  fis >> num_atoms;
  // A = 1.0 Angstrom (basic length-scale)
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  double basis_xx, basis_xy, basis_xz,
      basis_yx, basis_yy, basis_yz,
      basis_zx, basis_zy, basis_zz;
  // "H0(1,1) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_xx;
  // "H0(1,2) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_xy;
  // "H0(1,3) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_xz;
  // "H0(2,1) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_yx;
  // "H0(2,2) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_yy;
  // "H0(2,3) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_yz;
  // "H0(3,1) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_zx;
  // "H0(3,2) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_zy;
  // "H0(3,3) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_zz;
  // finish this line
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // .NO_VELOCITY.
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // "entry_count = 3"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto basis = Eigen::Matrix3d{{basis_xx, basis_xy, basis_xz},
                               {basis_yx, basis_yy, basis_yz},
                               {basis_zx, basis_zy, basis_zz}};
  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);

  double mass;
  std::string type;
  Eigen::Vector3d relative_position;

  std::vector<std::vector<size_t> > first_neighbors_adjacency_list,
      second_neighbors_adjacency_list, third_neighbors_adjacency_list;

  for (size_t id = 0; id < num_atoms; ++id) {
    fis >> mass >> type >> relative_position(0) >> relative_position(1) >> relative_position(2);
    atom_vector.emplace_back(type);
    relative_position_matrix.col(static_cast<int>(id)) = relative_position;
  }
  Config config_out = Config{basis, relative_position_matrix, atom_vector};
  config_out.ReassignLattice();
  return config_out;
}

Config Config::ReadPoscar(const std::string &filename) {
  std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
  if (!ifs) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  boost::iostreams::filtering_istream fis;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fis.push(boost::iostreams::gzip_decompressor());
  } else if (boost::filesystem::path(filename).extension() == ".bz2") {
    fis.push(boost::iostreams::bzip2_decompressor());
  }
  fis.push(ifs);

  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // #comment
  double scale;
  fis >> scale; // scale factor, usually which is 1.0
  Eigen::Matrix3d basis;
  fis >> basis(0, 0) >> basis(0, 1) >> basis(0, 2); // lattice vector a
  fis >> basis(1, 0) >> basis(1, 1) >> basis(1, 2); // lattice vector b
  fis >> basis(2, 0) >> basis(2, 1) >> basis(2, 2); // lattice vector c
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line
  basis *= scale;
  auto inverse_basis = basis.inverse();

  std::string buffer;
  getline(fis, buffer);
  std::istringstream element_iss(buffer);
  getline(fis, buffer);
  std::istringstream count_iss(buffer);

  std::string element;
  size_t num_atoms;
  std::vector<std::pair<std::string, size_t>> elements_counts;
  while (element_iss >> element && count_iss >> num_atoms) {
    elements_counts.emplace_back(element, num_atoms);
  }
  getline(fis, buffer);
  bool relative_option =
      buffer[0] != 'C' && buffer[0] != 'c' && buffer[0] != 'K' && buffer[0] != 'k';

  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);

  size_t id_count = 0;
  double position_X, position_Y, position_Z;
  for (const auto &[element_symbol, count] : elements_counts) {
    for (size_t j = 0; j < count; ++j) {
      fis >> position_X >> position_Y >> position_Z;
      fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      if (relative_option) {
        relative_position_matrix.col(static_cast<int>(id_count)) =
            Eigen::Vector3d(position_X, position_Y, position_Z);
      } else {
        relative_position_matrix.col(static_cast<int>(id_count)) =
            inverse_basis * Eigen::Vector3d(position_X, position_Y, position_Z);
      }
      atom_vector.emplace_back(element_symbol);
      ++id_count;
    }
  }
  Config config_out = Config{basis, relative_position_matrix, atom_vector};
  config_out.ReassignLattice();
  return config_out;
}

void Config::WriteConfig(const std::string &filename) const {
  WriteConfigExtended(filename, {});
}

void Config::WriteConfigExtended(
    const std::string &filename,
    const std::map<std::string, std::vector<double>> &auxiliary_lists) const {
  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fos.push(boost::iostreams::gzip_compressor());
  } else if (boost::filesystem::path(filename).extension() == ".bz2") {
    fos.push(boost::iostreams::bzip2_compressor());
  }
  fos.push(ofs);

  fos.precision(8);
  fos << "Number of particles = " << GetNumAtoms() << '\n';
  fos << "A = 1.0 Angstrom (basic length-scale)\n";
  fos << "H0(1,1) = " << basis_(0, 0) << " A\n";
  fos << "H0(1,2) = " << basis_(0, 1) << " A\n";
  fos << "H0(1,3) = " << basis_(0, 2) << " A\n";
  fos << "H0(2,1) = " << basis_(1, 0) << " A\n";
  fos << "H0(2,2) = " << basis_(1, 1) << " A\n";
  fos << "H0(2,3) = " << basis_(1, 2) << " A\n";
  fos << "H0(3,1) = " << basis_(2, 0) << " A\n";
  fos << "H0(3,2) = " << basis_(2, 1) << " A\n";
  fos << "H0(3,3) = " << basis_(2, 2) << " A\n";
  fos << ".NO_VELOCITY.\n";
  fos << "entry_count = " << 3 + auxiliary_lists.size() << "\n";
  size_t auxiliary_index = 0;
  for (const auto &auxiliary_list : auxiliary_lists) {
    fos << "auxiliary[" << auxiliary_index << "] = " << auxiliary_list.first << " [reduced unit]\n";
    ++auxiliary_index;
  }
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n", "", "", "", "");
  for (size_t it = 0; it < atom_vector_.size(); ++it) {
    const auto &atom = atom_vector_[it];
    const auto &relative_position
        = relative_position_matrix_.col(static_cast<int>(atom_to_lattice_hashmap_.at(it)));
    fos << atom.GetMass() << '\n'
        << atom.GetElementString() << '\n'
        << relative_position.transpose().format(fmt);
    for (const auto &[key, auxiliary_list] : auxiliary_lists) {
      fos << ' ' << auxiliary_list.at(it);
    }
    fos << std::endl;
  }
}

void Config::WriteXyzExtended(const std::string &filename,
                              const std::map<std::string, VectorVariant> &auxiliary_lists,
                              const std::map<std::string, ValueVariant> &global_list) const {
  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fos.push(boost::iostreams::gzip_compressor());
  } else if (boost::filesystem::path(filename).extension() == ".bz2") {
    fos.push(boost::iostreams::bzip2_compressor());
  }
  fos.push(ofs);
  fos.precision(8);
  fos << GetNumAtoms() << '\n';
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
  fos << "Lattice=\"" << basis_.format(fmt) << "\" ";
  fos << "pbc=\"" << periodic_boundary_condition_[0] << " " << periodic_boundary_condition_[1] << " "
      << periodic_boundary_condition_[2] << "\" ";
  for (const auto &[key, value] : global_list) {
    fos << key << "=";
    std::visit([&fos](const auto &val) { fos << val << ' '; }, value);
  }
  fos << "Properties=species:S:1:pos:R:3";
  for (const auto &[key, auxiliary_list] : auxiliary_lists) {
    fos << ":" << key << ":";
    std::any ret;
    std::visit([&ret](const auto &vec) { ret = vec.at(0); }, auxiliary_list);
    if (ret.type() == typeid(int) || ret.type() == typeid(size_t)) {
      fos << "I:1";
    } else if (ret.type() == typeid(double)) {
      fos << "R:1";
    } else if (ret.type() == typeid(std::string)) {
      fos << "S:1";
    } else if (ret.type() == typeid(Eigen::Vector3d) && std::any_cast<Eigen::Vector3d>(ret).size() == 3) {
      fos << "R:3";
    } else {
      throw std::runtime_error("Unsupported type");
    }
  }
  fos << '\n';

  for (size_t it = 0; it < atom_vector_.size(); ++it) {
    const auto &atom = atom_vector_[it];
    const auto &cartesian_position = cartesian_position_matrix_.col(static_cast<int>(atom_to_lattice_hashmap_.at(it)));
    fos << atom.GetElementString() << ' '
        << cartesian_position.transpose().format(fmt) << ' ';

    for (const auto &[key, auxiliary_list] : auxiliary_lists) {
      std::any ret;
      std::visit([&ret, it](const auto &vec) { ret = vec.at(it); }, auxiliary_list);
      if (ret.type() == typeid(int)) {
        fos << std::any_cast<int>(ret) << ' ';
      } else if (ret.type() == typeid(size_t)) {
        fos << std::any_cast<size_t>(ret) << ' ';
      } else if (ret.type() == typeid(double)) {
        fos << std::any_cast<double>(ret) << ' ';
      } else if (ret.type() == typeid(std::string)) {
        fos << std::any_cast<std::string>(ret) << ' ';
      } else if (ret.type() == typeid(Eigen::Vector3d)) {
        fos << std::any_cast<Eigen::Vector3d>(ret).format(fmt) << ' ';
      } else {
        throw std::runtime_error("Unsupported type");
      }
    }
    fos << std::endl;
  }
}
