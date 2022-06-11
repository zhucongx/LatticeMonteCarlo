#include "Config.h"

#include <utility>
#include <boost/functional/hash.hpp>

namespace cfg {
Config::Config() = default;
Config::Config(const Matrix_t &basis,
               std::vector<Lattice> lattice_vector,
               std::vector<Atom> atom_vector)
    : basis_(basis),
      lattice_vector_(std::move(lattice_vector)),
      atom_vector_(std::move(atom_vector)) {

  if (lattice_vector_.size() != atom_vector_.size()) {
    std::cerr << "Warning: lattice vector and atom vector size not match" << std::endl;
  }
  for (size_t i = 0; i < lattice_vector_.size(); ++i) {
    auto lattice_id = lattice_vector_.at(i).GetId();
    auto atom_id = atom_vector_.at(i).GetId();
    lattice_to_atom_hashmap_.emplace(lattice_id, atom_id);
    atom_to_lattice_hashmap_.emplace(atom_id, lattice_id);
  }
  UpdateNeighbors();
}
size_t Config::GetNumAtoms() const {
  return atom_vector_.size();
}
const Matrix_t &Config::GetBasis() const {
  return basis_;
}
const std::unordered_map<size_t, size_t> &Config::GetLatticeToAtomHashmap() const {
  return lattice_to_atom_hashmap_;
}
const std::unordered_map<size_t, size_t> &Config::GetAtomToLatticeHashmap() const {
  return atom_to_lattice_hashmap_;
}
const std::vector<Lattice> &Config::GetLatticeVector() const {
  return lattice_vector_;
}
const std::vector<Atom> &Config::GetAtomVector() const {
  return atom_vector_;
}
const std::vector<std::vector<size_t> > &Config::GetFirstNeighborsAdjacencyList() const {
  return first_neighbors_adjacency_list_;
}
const std::vector<std::vector<size_t> > &Config::GetSecondNeighborsAdjacencyList() const {
  return second_neighbors_adjacency_list_;
}
const std::vector<std::vector<size_t> > &Config::GetThirdNeighborsAdjacencyList() const {
  return third_neighbors_adjacency_list_;
}
std::vector<size_t> Config::GetFirstNeighborsAtomIdVectorOfAtom(size_t atom_id) const {
  auto lattice_id = atom_to_lattice_hashmap_.at(atom_id);
  std::vector<size_t> first_neighbors_atom_id_vector;
  first_neighbors_atom_id_vector.reserve(constants::kNumFirstNearestNeighbors);
  for (auto neighbor_lattice_id: first_neighbors_adjacency_list_[lattice_id]) {
    first_neighbors_atom_id_vector.push_back(lattice_to_atom_hashmap_.at(neighbor_lattice_id));
  }
  return first_neighbors_atom_id_vector;
}
size_t Config::GetAtomIdFromLatticeId(size_t lattice_id) const {
  return lattice_to_atom_hashmap_.at(lattice_id);
}
size_t Config::GetLatticeIdFromAtomId(size_t atom_id) const {
  return atom_to_lattice_hashmap_.at(atom_id);
}
Element Config::GetElementAtLatticeId(size_t lattice_id) const {
  auto atom_id = lattice_to_atom_hashmap_.at(lattice_id);
  return atom_vector_[atom_id].GetElement();
}
std::set<Element> Config::GetElementSetWithoutVacancy() const {
  std::set<Element> res;
  for (const auto &atom: atom_vector_) {
    if (atom.GetElement() == ElementName::X) { continue; }
    res.insert(atom.GetElement());
  }
  return res;
}
size_t Config::GetStateHash() const {
  size_t seed = 0;
  for (size_t i = 0; i < GetNumAtoms(); ++i) {
    boost::hash_combine(seed, lattice_to_atom_hashmap_.at(i));
  }
  return seed;
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
void Config::LatticeJump(const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  const auto [lattice_id_lhs, lattice_id_rhs] = lattice_id_jump_pair;
  const auto atom_id_lhs = lattice_to_atom_hashmap_.at(lattice_id_lhs);
  const auto atom_id_rhs = lattice_to_atom_hashmap_.at(lattice_id_rhs);

  atom_to_lattice_hashmap_.at(atom_id_lhs) = lattice_id_rhs;
  atom_to_lattice_hashmap_.at(atom_id_rhs) = lattice_id_lhs;
  lattice_to_atom_hashmap_.at(lattice_id_lhs) = atom_id_rhs;
  lattice_to_atom_hashmap_.at(lattice_id_rhs) = atom_id_lhs;
}
Config Config::ReadCfg(const std::string &filename) {
  std::ifstream ifs(filename, std::ifstream::in);
  // "Number of particles = %i"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  size_t num_atoms;
  ifs >> num_atoms;
  // A = 1.0 Angstrom (basic length-scale)
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  double basis_xx, basis_xy, basis_xz,
      basis_yx, basis_yy, basis_yz,
      basis_zx, basis_zy, basis_zz;
  // "H0(1,1) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_xx;
  // "H0(1,2) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_xy;
  // "H0(1,3) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_xz;
  // "H0(2,1) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_yx;
  // "H0(2,2) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_yy;
  // "H0(2,3) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_yz;
  // "H0(3,1) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_zx;
  // "H0(3,2) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_zy;
  // "H0(3,3) = %lf A"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  ifs >> basis_zz;
  // finish this line
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // .NO_VELOCITY.
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // "entry_count = 3"
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto basis = Matrix_t{{{basis_xx, basis_xy, basis_xz},
                         {basis_yx, basis_yy, basis_yz},
                         {basis_zx, basis_zy, basis_zz}}};
  std::vector<Atom> atom_vector;
  atom_vector.reserve(num_atoms);
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(num_atoms);

  double mass;
  std::string type;
  Vector_t relative_position;
  for (size_t id = 0; id < num_atoms; ++id) {
    ifs >> mass >> type >> relative_position;
    atom_vector.emplace_back(id, type);
    lattice_vector.emplace_back(id, relative_position * basis, relative_position);
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return Config{basis, lattice_vector, atom_vector};
}
void Config::WriteCfg(const std::string &filename, bool neighbors_info) const {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << "Number of particles = " << GetNumAtoms() << '\n';
  ofs << "A = 1.0 Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << basis_[kXDimension][kXDimension] << " A\n";
  ofs << "H0(1,2) = " << basis_[kXDimension][kYDimension] << " A\n";
  ofs << "H0(1,3) = " << basis_[kXDimension][kZDimension] << " A\n";
  ofs << "H0(2,1) = " << basis_[kYDimension][kXDimension] << " A\n";
  ofs << "H0(2,2) = " << basis_[kYDimension][kYDimension] << " A\n";
  ofs << "H0(2,3) = " << basis_[kYDimension][kZDimension] << " A\n";
  ofs << "H0(3,1) = " << basis_[kZDimension][kXDimension] << " A\n";
  ofs << "H0(3,2) = " << basis_[kZDimension][kYDimension] << " A\n";
  ofs << "H0(3,3) = " << basis_[kZDimension][kZDimension] << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom: atom_vector_) {
    size_t lattice_id = atom_to_lattice_hashmap_.at(atom.GetId());
    const auto &lattice = lattice_vector_[lattice_id];
    ofs << atom.GetMass() << '\n'
        << atom.GetElementString() << '\n'
        << lattice.GetRelativePosition();
    if (neighbors_info) {
      ofs << " # ";
      for (auto neighbor_lattice_index: first_neighbors_adjacency_list_[lattice_id]) {
        ofs << lattice_to_atom_hashmap_.at(neighbor_lattice_index) << ' ';
      }
      for (auto neighbor_lattice_index: second_neighbors_adjacency_list_[lattice_id]) {
        ofs << lattice_to_atom_hashmap_.at(neighbor_lattice_index) << ' ';
      }
      for (auto neighbor_lattice_index: third_neighbors_adjacency_list_[lattice_id]) {
        ofs << lattice_to_atom_hashmap_.at(neighbor_lattice_index) << ' ';
      }
    }
    ofs << '\n';
    ofs << std::flush;
  }
}
Config Config::ReadMap(const std::string &lattice_filename,
                       const std::string &element_filename,
                       const std::string &map_filename) {
  Config config;
  std::ifstream ifs_lattice(lattice_filename, std::ifstream::in);
  size_t num_atoms;
  ifs_lattice >> num_atoms;
  ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  Matrix_t basis;
  ifs_lattice >> basis;
  config.basis_ = basis;
  config.InitializeNeighborsList(num_atoms);

  config.lattice_vector_.reserve(num_atoms);
  Vector_t relative_position;
  size_t neighbor_id;
  for (size_t lattice_id = 0; lattice_id < num_atoms; ++lattice_id) {
    ifs_lattice >> relative_position;
    config.lattice_vector_.emplace_back(lattice_id, relative_position * basis, relative_position);
    ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '#');
    for (size_t i = 0; i < constants::kNumFirstNearestNeighbors; ++i) {
      ifs_lattice >> neighbor_id;
      config.first_neighbors_adjacency_list_[lattice_id].push_back(neighbor_id);
    }
    for (size_t i = 0; i < constants::kNumSecondNearestNeighbors; ++i) {
      ifs_lattice >> neighbor_id;
      config.second_neighbors_adjacency_list_[lattice_id].push_back(neighbor_id);
    }
    for (size_t i = 0; i < constants::kNumThirdNearestNeighbors; ++i) {
      ifs_lattice >> neighbor_id;
      config.third_neighbors_adjacency_list_[lattice_id].push_back(neighbor_id);
    }
    ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  std::ifstream ifs_element(element_filename, std::ifstream::in);
  std::string type;
  config.atom_vector_.reserve(num_atoms);
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id) {
    ifs_element >> type;
    config.atom_vector_.emplace_back(atom_id, type);
    ifs_element.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  std::ifstream ifs_map(map_filename, std::ifstream::in);
  size_t lattice_id;
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id) {
    ifs_map >> lattice_id;
    config.lattice_to_atom_hashmap_.emplace(lattice_id, atom_id);
    config.atom_to_lattice_hashmap_.emplace(atom_id, lattice_id);
    ifs_map.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  return config;
}
void Config::WriteLattice(const std::string &filename) const {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << GetNumAtoms() << " positions in total" << '\n';
  ofs << basis_ << std::endl;
  for (size_t i = 0; i < lattice_vector_.size(); ++i) {
    ofs << lattice_vector_[i].GetRelativePosition();
    ofs << " # ";
    for (auto neighbor_lattice_index: first_neighbors_adjacency_list_[i]) {
      ofs << lattice_to_atom_hashmap_.at(neighbor_lattice_index) << ' ';
    }
    for (auto neighbor_lattice_index: second_neighbors_adjacency_list_[i]) {
      ofs << lattice_to_atom_hashmap_.at(neighbor_lattice_index) << ' ';
    }
    for (auto neighbor_lattice_index: third_neighbors_adjacency_list_[i]) {
      ofs << lattice_to_atom_hashmap_.at(neighbor_lattice_index) << ' ';
    }
    ofs << std::endl;
  }
}
void Config::WriteElement(const std::string &filename) const {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  for (const auto &atom: atom_vector_) {
    ofs << atom.GetElementString() << std::endl;
  }
}
void Config::WriteMap(const std::string &filename) const {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  for (const auto &atom: atom_vector_) {
    size_t lattice_id = atom_to_lattice_hashmap_.at(atom.GetId());
    ofs << lattice_id << std::endl;
  }
}
void Config::ConvertRelativeToCartesian() {
  for (auto &lattice: lattice_vector_) {
    lattice.SetCartesianPosition(lattice.GetRelativePosition() * basis_);
  }
}
void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix(basis_);
  for (auto &lattice: lattice_vector_) {
    lattice.SetRelativePosition(lattice.GetCartesianPosition() * inverse_basis);
  }
}
void Config::InitializeNeighborsList(size_t num_atoms) {
  first_neighbors_adjacency_list_.resize(num_atoms);
  for (auto &neighbor_list: first_neighbors_adjacency_list_) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumFirstNearestNeighbors);
  }
  second_neighbors_adjacency_list_.resize(num_atoms);
  for (auto &neighbor_list: second_neighbors_adjacency_list_) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumSecondNearestNeighbors);
  }
  third_neighbors_adjacency_list_.resize(num_atoms);
  for (auto &neighbor_list: third_neighbors_adjacency_list_) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumThirdNearestNeighbors);
  }
}
void Config::UpdateNeighbors() {
  InitializeNeighborsList(GetNumAtoms());
  const double first_r_cutoff_square = std::pow(constants::kFirstNearestNeighborsCutoff, 2);
  const double second_r_cutoff_square = std::pow(constants::kSecondNearestNeighborsCutoff, 2);
  const double third_r_cutoff_square = std::pow(constants::kThirdNearestNeighborsCutoff, 2);
  for (auto it1 = atom_vector_.begin(); it1 != atom_vector_.end(); ++it1) {
    for (auto it2 = atom_vector_.begin(); it2 != it1; ++it2) {
      auto first_lattice_id = atom_to_lattice_hashmap_[it1->GetId()];
      auto second_lattice_id = atom_to_lattice_hashmap_[it2->GetId()];
      Vector_t absolute_distance_vector =
          GetRelativeDistanceVectorLattice(lattice_vector_[first_lattice_id],
                                           lattice_vector_[second_lattice_id]) * basis_;
      if (std::abs(absolute_distance_vector[kXDimension])
          > constants::kNearNeighborsCutoff) { continue; }
      if (std::abs(absolute_distance_vector[kYDimension])
          > constants::kNearNeighborsCutoff) { continue; }
      if (std::abs(absolute_distance_vector[kZDimension])
          > constants::kNearNeighborsCutoff) { continue; }
      const double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square < first_r_cutoff_square) {
        first_neighbors_adjacency_list_[first_lattice_id].push_back(second_lattice_id);
        first_neighbors_adjacency_list_[second_lattice_id].push_back(first_lattice_id);
      } else if (absolute_distance_square < second_r_cutoff_square) {
        second_neighbors_adjacency_list_[first_lattice_id].push_back(second_lattice_id);
        second_neighbors_adjacency_list_[second_lattice_id].push_back(first_lattice_id);
      } else if (absolute_distance_square < third_r_cutoff_square) {
        third_neighbors_adjacency_list_[first_lattice_id].push_back(second_lattice_id);
        third_neighbors_adjacency_list_[second_lattice_id].push_back(first_lattice_id);
      }
    }
  }
}

Vector_t GetLatticePairCenter(const Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  Vector_t center_position;
  for (const auto kDim: All_Dimensions) {
    double first_relative =
        config.GetLatticeVector()[lattice_id_jump_pair.first].GetRelativePosition()[kDim];
    const double second_relative =
        config.GetLatticeVector()[lattice_id_jump_pair.second].GetRelativePosition()[kDim];

    double distance = first_relative - second_relative;
    int period = static_cast<int>(distance / 0.5);
    // make sure distance is the range (0, 0.5)
    while (period != 0) {
      first_relative -= static_cast<double>(period);
      distance = first_relative - second_relative;
      period = static_cast<int>(distance / 0.5);
    }
    center_position[kDim] = 0.5 * (first_relative + second_relative);
  }
  return center_position;
}
Matrix_t GetLatticePairRotationMatrix(const Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  const auto &first_lattice = config.GetLatticeVector()[lattice_id_jump_pair.first];
  const auto &second_lattice = config.GetLatticeVector()[lattice_id_jump_pair.second];

  const Vector_t
      pair_direction = Normalize(GetRelativeDistanceVectorLattice(first_lattice, second_lattice));
  Vector_t vertical_vector{};
  for (const auto index: config.GetFirstNeighborsAdjacencyList().at(lattice_id_jump_pair.first)) {
    const Vector_t jump_vector =
        GetRelativeDistanceVectorLattice(first_lattice, config.GetLatticeVector()[index]);
    const double dot_prod = Dot(pair_direction, jump_vector);
    if (std::abs(dot_prod) < 1e-6) {
      vertical_vector = Normalize(jump_vector);
      break;
    }
  }
  // The third row is normalized since it is a cross product of two normalized vectors.
  // We use transposed matrix here because transpose of an orthogonal matrix equals its inverse
  return TransposeMatrix({pair_direction, vertical_vector,
                          Cross(pair_direction, vertical_vector)});
}
void RotateLatticeVector(std::vector<Lattice> &lattice_list,
                         const Matrix_t &rotation_matrix) {
  const auto move_distance_after_rotation = Vector_t{0.5, 0.5, 0.5}
      - (Vector_t{0.5, 0.5, 0.5} * rotation_matrix);
  for (auto &lattice: lattice_list) {
    auto relative_position = lattice.GetRelativePosition();
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= ElementFloor(relative_position);
    lattice.SetRelativePosition(relative_position);
  }
}
size_t GetVacancyAtomIndex(const Config &config) {
  return config.GetAtomIdFromLatticeId(GetVacancyLatticeIndex(config));
}

size_t GetVacancyLatticeIndex(const Config &config) {
  const auto &atom_vector = config.GetAtomVector();
  auto it = std::find_if(atom_vector.cbegin(),
                         atom_vector.cend(),
                         [](const auto &atom) {
                           return atom.GetElement() == ElementName::X;
                         });
  if (it != atom_vector.end()) {
    return config.GetLatticeIdFromAtomId(it->GetId());
  } else {
    std::cerr << "Warning: vacancy not found" << std::endl;
  }
  return 0;
}

std::unordered_set<size_t> GetNeighborsLatticeIdSetOfLattice(
    const Config &config, size_t lattice_id) {
  std::unordered_set<size_t> near_neighbors_hashset;
  std::copy(config.GetFirstNeighborsAdjacencyList().at(lattice_id).begin(),
            config.GetFirstNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset,
                          near_neighbors_hashset.begin()));
  std::copy(config.GetSecondNeighborsAdjacencyList().at(lattice_id).begin(),
            config.GetSecondNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset,
                          near_neighbors_hashset.begin()));
  std::copy(config.GetThirdNeighborsAdjacencyList().at(lattice_id).begin(),
            config.GetThirdNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset,
                          near_neighbors_hashset.begin()));
  return near_neighbors_hashset;
}
std::unordered_set<size_t> GetNeighborsLatticeIdSetOfJumpPair(
    const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {

  std::unordered_set<size_t> near_neighbors_hashset;
  for (const auto lattice_id: {lattice_id_jump_pair.first, lattice_id_jump_pair.second}) {
    std::copy(config.GetFirstNeighborsAdjacencyList().at(lattice_id).begin(),
              config.GetFirstNeighborsAdjacencyList().at(lattice_id).end(),
              std::inserter(near_neighbors_hashset,
                            near_neighbors_hashset.begin()));
    std::copy(config.GetSecondNeighborsAdjacencyList().at(lattice_id).begin(),
              config.GetSecondNeighborsAdjacencyList().at(lattice_id).end(),
              std::inserter(near_neighbors_hashset,
                            near_neighbors_hashset.begin()));
    std::copy(config.GetThirdNeighborsAdjacencyList().at(lattice_id).begin(),
              config.GetThirdNeighborsAdjacencyList().at(lattice_id).end(),
              std::inserter(near_neighbors_hashset,
                            near_neighbors_hashset.begin()));
  }
  return near_neighbors_hashset;
}
int FindDistanceLabelBetweenLattice(size_t index1, size_t index2, const Config &config) {
  const auto &first_neighbors_adjacency_list = config.GetFirstNeighborsAdjacencyList()[index1];
  const auto &second_neighbors_adjacency_list = config.GetSecondNeighborsAdjacencyList()[index1];
  const auto &third_neighbors_adjacency_list = config.GetThirdNeighborsAdjacencyList()[index1];
  if (std::find(first_neighbors_adjacency_list.begin(),
                first_neighbors_adjacency_list.end(),
                index2) != first_neighbors_adjacency_list.end()) {
    return 1;
  }
  if (std::find(second_neighbors_adjacency_list.begin(),
                second_neighbors_adjacency_list.end(),
                index2) != second_neighbors_adjacency_list.end()) {
    return 2;
  }
  if (std::find(third_neighbors_adjacency_list.begin(),
                third_neighbors_adjacency_list.end(),
                index2) != third_neighbors_adjacency_list.end()) {
    return 3;
  }
  return -1;

}

} // namespace cfg
