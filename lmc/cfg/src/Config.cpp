#include "Config.h"

#include <algorithm>
#include <any>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <chrono>
#include <omp.h>
#include <random>
#include <utility>

namespace cfg {
Config::Config() = default;

Config::Config(const Matrix_t &basis,
               std::vector<Lattice> lattice_vector,
               std::vector<Atom> atom_vector,
               bool update_neighbor)
    : basis_(basis),
      lattice_vector_(std::move(lattice_vector)),
      atom_vector_(std::move(atom_vector)) {
  map_shift_list_.resize(atom_vector_.size(), {0, 0, 0});
  if (lattice_vector_.size() != atom_vector_.size()) {
    throw std::runtime_error("Lattice vector and atom vector size do not match");
  }
  for (size_t i = 0; i < lattice_vector_.size(); ++i) {
    auto lattice_id = lattice_vector_.at(i).GetId();
    auto atom_id = atom_vector_.at(i).GetId();
    lattice_to_atom_hashmap_.emplace(lattice_id, atom_id);
    atom_to_lattice_hashmap_.emplace(atom_id, lattice_id);
  }
  if (update_neighbor) {
    UpdateNeighbors();
  }
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

const std::vector<std::array<int, 3>> &Config::GetMapShiftList() const {
  return map_shift_list_;
}

const std::vector<std::vector<size_t>> &Config::GetFirstNeighborsAdjacencyList() const {
  return first_neighbors_adjacency_list_;
}

const std::vector<std::vector<size_t>> &Config::GetSecondNeighborsAdjacencyList() const {
  return second_neighbors_adjacency_list_;
}

const std::vector<std::vector<size_t>> &Config::GetThirdNeighborsAdjacencyList() const {
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

std::vector<size_t> Config::GetSecondNeighborsAtomIdVectorOfAtom(size_t atom_id) const {
  const auto lattice_id = atom_to_lattice_hashmap_.at(atom_id);
  std::vector<size_t> second_neighbors_atom_id_vector;
  second_neighbors_atom_id_vector.reserve(constants::kNumSecondNearestNeighbors);
  for (auto neighbor_lattice_id: second_neighbors_adjacency_list_[lattice_id]) {
    second_neighbors_atom_id_vector.push_back(lattice_to_atom_hashmap_.at(neighbor_lattice_id));
  }
  return second_neighbors_atom_id_vector;
}

std::vector<size_t> Config::GetThirdNeighborsAtomIdVectorOfAtom(size_t atom_id) const {
  const auto lattice_id = atom_to_lattice_hashmap_.at(atom_id);
  std::vector<size_t> third_neighbors_atom_id_vector;
  third_neighbors_atom_id_vector.reserve(constants::kNumThirdNearestNeighbors);
  for (auto neighbor_lattice_id: third_neighbors_adjacency_list_[lattice_id]) {
    third_neighbors_atom_id_vector.push_back(lattice_to_atom_hashmap_.at(neighbor_lattice_id));
  }
  return third_neighbors_atom_id_vector;
}

size_t Config::GetAtomIdFromLatticeId(size_t lattice_id) const {
  return lattice_to_atom_hashmap_.at(lattice_id);
}

size_t Config::GetLatticeIdFromAtomId(size_t atom_id) const {
  return atom_to_lattice_hashmap_.at(atom_id);
}

Element Config::GetElementAtAtomId(size_t atom_id) const {
  return atom_vector_[atom_id].GetElement();
}

Element Config::GetElementAtLatticeId(size_t lattice_id) const {
  auto atom_id = lattice_to_atom_hashmap_.at(lattice_id);
  return atom_vector_[atom_id].GetElement();
}

std::set<Element> Config::GetElementSetWithoutVacancy() const {
  std::set<Element> res;
  for (const auto &atom: atom_vector_) {
    if (atom.GetElement() == ElementName::X) {
      continue;
    }
    res.insert(atom.GetElement());
  }
  return res;
}

double Config::GetTotalSoluteMass() const {
  const auto solvent_element = GetSolventElement();
  double total_mass = 0.0;

  for (const auto& atom : atom_vector_) {
    // Skip solvent atoms
    if (atom.GetElement() != solvent_element) {
      total_mass += atom.GetElement().GetMass();
    }
  }

  return total_mass;
}

Vector_t Config::GetSoluteCenterOfMass() const {
  const auto solvent_element = GetSolventElement();
  Vector_t mass_center{};
  Vector_t sum_cos_theta{};
  Vector_t sum_sin_theta{};
  double sum_mass = 0.0;

  for (const auto& atom : atom_vector_) {
    // Skip solvent atoms and vacancy
    if (atom.GetElement() == solvent_element || atom.GetElement() == ElementName::X) {
      continue;
    }

    const auto& lattice_id = atom_to_lattice_hashmap_.at(atom.GetId());
    const auto& relative_position = lattice_vector_[lattice_id].GetRelativePosition();
    const double mass = atom.GetElement().GetMass();

    sum_mass += mass;
    for (const auto kDim : All_Dimensions) {
      const double theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta) * mass;
      sum_sin_theta[kDim] += std::sin(theta) * mass;
    }
  }

  // If there are no solute atoms, return zero vector
  if (std::abs(sum_mass) < kEpsilon) {
    return mass_center;
  }

  const auto cos_theta_bar = sum_cos_theta / sum_mass;
  const auto sin_theta_bar = sum_sin_theta / sum_mass;

  for (const auto kDim : All_Dimensions) {
    const double theta_bar = std::atan2(-sin_theta_bar[kDim], -cos_theta_bar[kDim]) + M_PI;
    mass_center[kDim] = theta_bar / (2 * M_PI);
  }

  // Convert from relative to Cartesian coordinates
  return mass_center * basis_;
}

std::map<Element, std::vector<size_t>> Config::GetElementAtomIdVectorMap() const {
  std::map<Element, std::vector<size_t>> element_list_map;
  for (const auto &atom: atom_vector_) {
    element_list_map[atom.GetElement()].push_back(atom.GetId());
  }
  return element_list_map;
}

std::map<Element, size_t> Config::GetElementCountMap() const {
  std::map<Element, size_t> element_count_map;
  for (const auto &atom: atom_vector_) {
    element_count_map[atom.GetElement()]++;
  }
  return element_count_map;
}

size_t Config::GetStateHash() const {
  size_t seed = 0;
  for (size_t i = 0; i < GetNumAtoms(); ++i) {
    boost::hash_combine(seed, lattice_to_atom_hashmap_.at(i));
  }
  return seed;
}

Vector_t Config::GetLatticePairCenter(const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  Vector_t center_position;
  for (const auto kDim: All_Dimensions) {
    double first_relative = GetLatticeVector()[lattice_id_jump_pair.first].GetRelativePosition()[kDim];
    const double second_relative = GetLatticeVector()[lattice_id_jump_pair.second].GetRelativePosition()[kDim];

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

Matrix_t Config::GetLatticePairRotationMatrix(const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto &first_lattice = GetLatticeVector()[lattice_id_jump_pair.first];
  const auto &second_lattice = GetLatticeVector()[lattice_id_jump_pair.second];

  const Vector_t pair_direction = Normalize(GetRelativeDistanceVectorLattice(first_lattice, second_lattice)*GetBasis());
  Vector_t vertical_vector{};
  for (const auto index: GetFirstNeighborsAdjacencyList().at(lattice_id_jump_pair.first)) {
    const Vector_t jump_vector = Normalize(GetRelativeDistanceVectorLattice(first_lattice, GetLatticeVector()[index])*GetBasis());
    const double dot_prod = Dot(pair_direction, jump_vector);
    if (std::abs(dot_prod) < kEpsilon) {
      vertical_vector = jump_vector;
      break;
    }
  }
  // The third row is normalized since it is a cross product of two normalized vectors.
  // We use transposed matrix here because transpose of an orthogonal matrix equals its inverse
  return TransposeMatrix({pair_direction, vertical_vector, Cross(pair_direction, vertical_vector)});
}

size_t Config::GetVacancyAtomId() const {
  return GetAtomIdFromLatticeId(GetVacancyLatticeId());
}

size_t Config::GetVacancyLatticeId() const {
  const auto &atom_vector = GetAtomVector();
  auto it = std::find_if(atom_vector.cbegin(), atom_vector.cend(), [](const auto &atom) {
    return atom.GetElement() == ElementName::X;
  });
  if (it != atom_vector.end()) {
    return GetLatticeIdFromAtomId(it->GetId());
  } else {
    throw std::runtime_error("vacancy not found");
    // std::cerr << "vacancy not found" << std::endl;
    // return GetLatticeIdFromAtomId(0);
  }
}

std::unordered_set<size_t> Config::GetNeighborsLatticeIdSetOfSite(size_t lattice_id) const {
  std::unordered_set<size_t> near_neighbors_hashset;
  near_neighbors_hashset.insert(lattice_id);
  std::copy(GetFirstNeighborsAdjacencyList().at(lattice_id).begin(),
            GetFirstNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
  std::copy(GetSecondNeighborsAdjacencyList().at(lattice_id).begin(),
            GetSecondNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
  std::copy(GetThirdNeighborsAdjacencyList().at(lattice_id).begin(),
            GetThirdNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
  return near_neighbors_hashset;
}
std::unordered_set<size_t> Config::GetNeighborsAtomIdSetOfAtom(size_t atom_id) const {
  std::unordered_set<size_t> near_neighbors_hashset;
  near_neighbors_hashset.insert(atom_id);
  const auto lattice_id = atom_to_lattice_hashmap_.at(atom_id);
  const auto lattice_id_set = GetNeighborsLatticeIdSetOfSite(lattice_id);
  for (const auto neighbor_lattice_id: lattice_id_set) {
    near_neighbors_hashset.insert(lattice_to_atom_hashmap_.at(neighbor_lattice_id));
  }
  return near_neighbors_hashset;
}
std::unordered_set<size_t>
Config::GetNeighborsLatticeIdSetOfPair(const std::pair<size_t, size_t> &lattice_id_pair) const {
  std::unordered_set<size_t> near_neighbors_hashset;
  for (const auto lattice_id: {lattice_id_pair.first, lattice_id_pair.second}) {
    near_neighbors_hashset.insert(lattice_id);
    std::copy(GetFirstNeighborsAdjacencyList().at(lattice_id).begin(),
              GetFirstNeighborsAdjacencyList().at(lattice_id).end(),
              std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
    std::copy(GetSecondNeighborsAdjacencyList().at(lattice_id).begin(),
              GetSecondNeighborsAdjacencyList().at(lattice_id).end(),
              std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
    std::copy(GetThirdNeighborsAdjacencyList().at(lattice_id).begin(),
              GetThirdNeighborsAdjacencyList().at(lattice_id).end(),
              std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
  }
  return near_neighbors_hashset;
}

int Config::FindDistanceLabelBetweenLattice(size_t lattice_id1, size_t lattice_id2) const {
  const auto &first_neighbors_adjacency_list = GetFirstNeighborsAdjacencyList()[lattice_id1];
  const auto &second_neighbors_adjacency_list = GetSecondNeighborsAdjacencyList()[lattice_id1];
  const auto &third_neighbors_adjacency_list = GetThirdNeighborsAdjacencyList()[lattice_id1];
  if (std::find(first_neighbors_adjacency_list.begin(), first_neighbors_adjacency_list.end(), lattice_id2) !=
      first_neighbors_adjacency_list.end()) {
    return 1;
  }
  if (std::find(second_neighbors_adjacency_list.begin(), second_neighbors_adjacency_list.end(), lattice_id2) !=
      second_neighbors_adjacency_list.end()) {
    return 2;
  }
  if (std::find(third_neighbors_adjacency_list.begin(), third_neighbors_adjacency_list.end(), lattice_id2) !=
      third_neighbors_adjacency_list.end()) {
    return 3;
  }
  return -1;
}

double Config::GetVacancyConcentration() const {
  size_t vacancy_count = 0;
  for (const auto &atom: GetAtomVector()) {
    if (atom.GetElement() == ElementName::X) {
      vacancy_count++;
    }
  }
  return static_cast<double>(vacancy_count) / static_cast<double>(GetNumAtoms());
}

double Config::GetSoluteConcentration(Element solvent_element) const {
  size_t solute_count = 0;
  for (const auto &atom: GetAtomVector()) {
    if (atom.GetElement() != solvent_element && atom.GetElement() != ElementName::X) {
      solute_count++;
    }
  }
  return static_cast<double>(solute_count) / static_cast<double>(GetNumAtoms());
}

std::map<Element, size_t> Config::GetLocalInfoOfLatticeId(size_t lattice_id, size_t shell_number) const {
  std::map<Element, size_t> element_count_map;
  for (const auto &element: GetElementSetWithoutVacancy()) {
    element_count_map.emplace(element, 0);
  }
  const std::vector<size_t> *neighbor_list = nullptr;
  switch (shell_number) {
    case 1: neighbor_list = &first_neighbors_adjacency_list_.at(lattice_id); break;
    case 2: neighbor_list = &second_neighbors_adjacency_list_.at(lattice_id); break;
    case 3: neighbor_list = &third_neighbors_adjacency_list_.at(lattice_id); break;
    default: throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
  }

  for (const auto neighbor_lattice_id: *neighbor_list) {
    element_count_map[GetElementAtLatticeId(neighbor_lattice_id)]++;
  }
  return element_count_map;
}
Element Config::GetSolventElement() const {
 auto element_atom_id_vector_map = GetElementAtomIdVectorMap();
 auto solvent_element = Element("X");
 size_t max_count = 0;
 for (const auto &[element, atom_id_vector]: element_atom_id_vector_map) {
   if (atom_id_vector.size() > max_count) {
     max_count = atom_id_vector.size();
     solvent_element = element;
   }
 }
 return solvent_element;
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

  const Vector_t displacement =
      lattice_vector_[lattice_id_rhs].GetRelativePosition() - lattice_vector_[lattice_id_lhs].GetRelativePosition();

  const auto image_change = ElementFloor(displacement + 0.5);

  for (const auto kDim : All_Dimensions) {
    map_shift_list_[atom_id_lhs][kDim] += static_cast<int>(image_change[kDim]);
    map_shift_list_[atom_id_rhs][kDim] -= static_cast<int>(image_change[kDim]);
  }

  atom_to_lattice_hashmap_.at(atom_id_lhs) = lattice_id_rhs;
  atom_to_lattice_hashmap_.at(atom_id_rhs) = lattice_id_lhs;
  lattice_to_atom_hashmap_.at(lattice_id_lhs) = atom_id_rhs;
  lattice_to_atom_hashmap_.at(lattice_id_rhs) = atom_id_lhs;
}

void Config::SetAtomElementTypeAtAtom(size_t atom_id, Element element) {
  atom_vector_.at(atom_id).SetElement(element);
}

void Config::SetAtomElementTypeAtLattice(size_t lattice_id, Element element) {
  atom_vector_.at(lattice_to_atom_hashmap_.at(lattice_id)).SetElement(element);
}

void Config::ReassignLatticeVector() {
  auto new_lattice_vector(lattice_vector_);
  std::sort(new_lattice_vector.begin(), new_lattice_vector.end(), [](const auto &lhs, const auto &rhs) -> bool {
    return lhs.GetRelativePosition() < rhs.GetRelativePosition();
  });
  std::unordered_map<size_t, size_t> old_lattice_id_to_new;
  for (size_t lattice_id = 0; lattice_id < GetNumAtoms(); ++lattice_id) {
    old_lattice_id_to_new.emplace(new_lattice_vector.at(lattice_id).GetId(), lattice_id);
    new_lattice_vector.at(lattice_id).SetId(lattice_id);
  }
  std::unordered_map<size_t, size_t> new_lattice_to_atom_hashmap, new_atom_to_lattice_hashmap;
  for (size_t atom_id = 0; atom_id < GetNumAtoms(); ++atom_id) {
    auto old_lattice_id = atom_to_lattice_hashmap_.at(atom_id);
    auto new_lattice_id = old_lattice_id_to_new.at(old_lattice_id);
    new_lattice_to_atom_hashmap.emplace(new_lattice_id, atom_id);
    new_atom_to_lattice_hashmap.emplace(atom_id, new_lattice_id);
  }

  std::vector<std::vector<size_t>> new_first_neighbors_adjacency_list, new_second_neighbors_adjacency_list,
      new_third_neighbors_adjacency_list;
  new_first_neighbors_adjacency_list.resize(GetNumAtoms());
  for (auto &neighbor_list: new_first_neighbors_adjacency_list) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumFirstNearestNeighbors);
  }
  new_second_neighbors_adjacency_list.resize(GetNumAtoms());
  for (auto &neighbor_list: new_second_neighbors_adjacency_list) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumSecondNearestNeighbors);
  }
  new_third_neighbors_adjacency_list.resize(GetNumAtoms());
  for (auto &neighbor_list: new_third_neighbors_adjacency_list) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumThirdNearestNeighbors);
  }

  for (size_t old_lattice_id = 0; old_lattice_id < GetNumAtoms(); ++old_lattice_id) {
    const auto &old_neighbor_id_vector = first_neighbors_adjacency_list_.at(old_lattice_id);
    const auto new_lattice_id = old_lattice_id_to_new.at(old_lattice_id);
    for (const auto old_neighbor_id: old_neighbor_id_vector) {
      new_first_neighbors_adjacency_list.at(new_lattice_id).push_back(old_lattice_id_to_new.at(old_neighbor_id));
    }
  }
  for (size_t old_lattice_id = 0; old_lattice_id < GetNumAtoms(); ++old_lattice_id) {
    const auto &old_neighbor_id_vector = second_neighbors_adjacency_list_.at(old_lattice_id);
    const auto new_lattice_id = old_lattice_id_to_new.at(old_lattice_id);
    for (const auto old_neighbor_id: old_neighbor_id_vector) {
      new_second_neighbors_adjacency_list.at(new_lattice_id).push_back(old_lattice_id_to_new.at(old_neighbor_id));
    }
  }
  for (size_t old_lattice_id = 0; old_lattice_id < GetNumAtoms(); ++old_lattice_id) {
    const auto &old_neighbor_id_vector = third_neighbors_adjacency_list_.at(old_lattice_id);
    const auto new_lattice_id = old_lattice_id_to_new.at(old_lattice_id);
    for (const auto old_neighbor_id: old_neighbor_id_vector) {
      new_third_neighbors_adjacency_list.at(new_lattice_id).push_back(old_lattice_id_to_new.at(old_neighbor_id));
    }
  }

  for (size_t lattice_id = 0; lattice_id < GetNumAtoms(); ++lattice_id) {
    std::sort(new_first_neighbors_adjacency_list.at(lattice_id).begin(),
              new_first_neighbors_adjacency_list.at(lattice_id).end());
    std::sort(new_second_neighbors_adjacency_list.at(lattice_id).begin(),
              new_second_neighbors_adjacency_list.at(lattice_id).end());
    std::sort(new_third_neighbors_adjacency_list.at(lattice_id).begin(),
              new_third_neighbors_adjacency_list.at(lattice_id).end());
  }
  lattice_vector_ = new_lattice_vector;
  lattice_to_atom_hashmap_ = new_lattice_to_atom_hashmap;
  atom_to_lattice_hashmap_ = new_atom_to_lattice_hashmap;
  first_neighbors_adjacency_list_ = new_first_neighbors_adjacency_list;
  second_neighbors_adjacency_list_ = new_second_neighbors_adjacency_list;
  third_neighbors_adjacency_list_ = new_third_neighbors_adjacency_list;
}

Config Config::ReadConfig(const std::string &filename) {
  std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
  if (!ifs) {
    throw std::runtime_error("Cannot open " + filename);
  }
  boost::iostreams::filtering_istream fis;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fis.push(boost::iostreams::gzip_decompressor());
  }
  fis.push(ifs);
  // "Number of particles = %i"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  size_t num_atoms;
  fis >> num_atoms;
  // A = 1.0 Angstrom (basic length-scale)
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  double basis_xx, basis_xy, basis_xz, basis_yx, basis_yy, basis_yz, basis_zx, basis_zy, basis_zz;
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
  if (fis.peek() == '.') {
    fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  // "entry_count = 3"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto basis =
      Matrix_t{{{basis_xx, basis_xy, basis_xz}, {basis_yx, basis_yy, basis_yz}, {basis_zx, basis_zy, basis_zz}}};
  std::vector<Atom> atom_vector;
  atom_vector.reserve(num_atoms);
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(num_atoms);

  double mass;
  std::string type;
  Vector_t relative_position;

  std::vector<std::vector<size_t>> first_neighbors_adjacency_list, second_neighbors_adjacency_list,
      third_neighbors_adjacency_list;

  for (size_t lattice_id = 0; lattice_id < num_atoms; ++lattice_id) {
    fis >> mass >> type >> relative_position;
    atom_vector.emplace_back(lattice_id, type);
    lattice_vector.emplace_back(lattice_id, relative_position * basis, relative_position);

    fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  return Config{basis, lattice_vector, atom_vector, true};
}

void Config::WriteConfig(const std::string &filename) const {
  WriteExtendedConfig(filename, {});
}

void Config::WriteExtendedConfig(const std::string &filename,
                                 const std::map<std::string, std::vector<double>> &auxiliary_lists) const {
  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fos.push(boost::iostreams::gzip_compressor());
  }
  fos.push(ofs);
  fos.precision(16);
  fos << "Number of particles = " << GetNumAtoms() << '\n';
  fos << "A = 1.0 Angstrom (basic length-scale)\n";
  fos << "H0(1,1) = " << basis_[0][0] << " A\n";
  fos << "H0(1,2) = " << basis_[0][1] << " A\n";
  fos << "H0(1,3) = " << basis_[0][2] << " A\n";
  fos << "H0(2,1) = " << basis_[1][0] << " A\n";
  fos << "H0(2,2) = " << basis_[1][1] << " A\n";
  fos << "H0(2,3) = " << basis_[1][2] << " A\n";
  fos << "H0(3,1) = " << basis_[2][0] << " A\n";
  fos << "H0(3,2) = " << basis_[2][1] << " A\n";
  fos << "H0(3,3) = " << basis_[2][2] << " A\n";
  fos << ".NO_VELOCITY.\n";
  fos << "entry_count = " << 3 + auxiliary_lists.size() << "\n";

  size_t auxiliary_index = 0;
  for (const auto &auxiliary_list: auxiliary_lists) {
    fos << "auxiliary[" << auxiliary_index << "] = " << auxiliary_list.first << " [reduced unit]\n";
    ++auxiliary_index;
  }

  for (size_t it = 0; it < atom_vector_.size(); ++it) {
    const auto &atom = atom_vector_[it];
    fos << atom.GetMass() << '\n'
        << atom.GetElementString() << '\n'
        << lattice_vector_[atom_to_lattice_hashmap_.at(atom.GetId())].GetRelativePosition();
    for (const auto &auxiliary_list: auxiliary_lists) {
      ofs << ' ' << auxiliary_list.second[it];
    }
    fos << std::endl;
  }
}

void Config::WriteExtendedXyz(const std::string &filename,
                              const std::map<std::string, VectorVariant> &auxiliary_lists,
                              const std::map<std::string, ValueVariant> &global_list) const {
  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  if (boost::filesystem::path(filename).extension() == ".gz") {
    fos.push(boost::iostreams::gzip_compressor());
  }
  fos.push(ofs);
  fos.precision(16);
  fos << GetNumAtoms() << '\n';
  fos << "Lattice=\"" << basis_[0] << " " << basis_[1] << " " << basis_[2] << "\" ";
  fos << "pbc=\"T T T\" ";
  for (const auto &[key, value]: global_list) {
    fos << key << "=";
    std::visit(
        [&fos](const auto &val) {
          using ValueType = std::decay_t<decltype(val)>;
          if constexpr (std::is_same_v<ValueType, std::vector<double>>) {
            fos << "\"";
            for (size_t i = 0; i < val.size(); ++i) {
              fos << val[i];
              if (i < val.size() - 1) {
                fos << " ";
              }
            }
            fos << "\" ";
          } else {
            fos << val << ' ';
          }
        },
        value);
  }
  fos << "Properties=species:S:1:pos:R:3";
  for (const auto &[key, auxiliary_list]: auxiliary_lists) {
    fos << ":" << key << ":";
    std::any ret;
    std::visit(
        [&ret](const auto &vec) {
          ret = vec.at(0);
        },
        auxiliary_list);
    if (ret.type() == typeid(int) || ret.type() == typeid(size_t)) {
      fos << "I:1";
    } else if (ret.type() == typeid(double)) {
      fos << "R:1";
    } else if (ret.type() == typeid(std::string)) {
      fos << "S:1";
    } else if (ret.type() == typeid(Vector_t)) {
      fos << "R:3";
    } else if (ret.type() == typeid(std::array<int, 3>)) {
      fos << "I:3";
    } else if (ret.type() == typeid(std::vector<double>)) {
      fos << "R:" << std::any_cast<std::vector<double>>(ret).size();
    } else if (ret.type() == typeid(std::vector<size_t>)) {
      fos << "I:" << std::any_cast<std::vector<size_t>>(ret).size();
    } else {
      throw std::runtime_error("Unsupported type");
    }
  }
  fos << '\n';

  for (size_t it = 0; it < GetNumAtoms(); ++it) {
    const auto &cartesian_position = lattice_vector_[atom_to_lattice_hashmap_.at(it)].GetCartesianPosition();
    fos << atom_vector_[it].GetElementString() << ' ' << cartesian_position << ' ';

    for (const auto &[key, auxiliary_list]: auxiliary_lists) {
      std::any ret;
      std::visit(
          [&ret, it](const auto &vec) {
            ret = vec.at(it);
          },
          auxiliary_list);
      if (ret.type() == typeid(int)) {
        fos << std::any_cast<int>(ret) << ' ';
      } else if (ret.type() == typeid(size_t)) {
        fos << std::any_cast<size_t>(ret) << ' ';
      } else if (ret.type() == typeid(double)) {
        fos << std::any_cast<double>(ret) << ' ';
      } else if (ret.type() == typeid(std::string)) {
        fos << std::any_cast<std::string>(ret) << ' ';
      } else if (ret.type() == typeid(Vector_t)) {
        fos << std::any_cast<Vector_t>(ret) << ' ';
      } else if (ret.type() == typeid(std::array<int, 3>)) {
        const auto &val = std::any_cast<std::array<int, 3>>(ret);
        fos << val[0] << ' ' << val[1] << ' ' << val[2] << ' ';
      } else if (ret.type() == typeid(std::vector<double>)) {
        for (const auto &val: std::any_cast<std::vector<double>>(ret)) {
          fos << val << ' ';
        }
      } else if (ret.type() == typeid(std::vector<size_t>)) {
        for (const auto &val: std::any_cast<std::vector<size_t>>(ret)) {
          fos << val << ' ';
        }
      } else {
        throw std::runtime_error("Unsupported type");
      }
    }
    fos << std::endl;
  }
}

Config Config::ReadMap(const std::string &lattice_filename,
                       const std::string &element_filename,
                       const std::string &map_filename) {
  Config config;
  std::ifstream ifs_lattice(lattice_filename, std::ifstream::in);
  if (!ifs_lattice.is_open()) {
    throw std::runtime_error("Cannot open " + lattice_filename);
  }
  size_t num_atoms;
  ifs_lattice >> num_atoms;
  ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  Matrix_t basis;
  ifs_lattice >> basis;
  config.basis_ = basis;
  config.InitializeNeighborsList(num_atoms);

  config.lattice_vector_.reserve(num_atoms);
  Vector_t relative_position;
  size_t neighbor_lattice_id;
  for (size_t lattice_id = 0; lattice_id < num_atoms; ++lattice_id) {
    ifs_lattice >> relative_position;
    config.lattice_vector_.emplace_back(lattice_id, relative_position * basis, relative_position);
    ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '#');
    for (size_t i = 0; i < constants::kNumFirstNearestNeighbors; ++i) {
      ifs_lattice >> neighbor_lattice_id;
      config.first_neighbors_adjacency_list_[lattice_id].push_back(neighbor_lattice_id);
    }
    for (size_t i = 0; i < constants::kNumSecondNearestNeighbors; ++i) {
      ifs_lattice >> neighbor_lattice_id;
      config.second_neighbors_adjacency_list_[lattice_id].push_back(neighbor_lattice_id);
    }
    for (size_t i = 0; i < constants::kNumThirdNearestNeighbors; ++i) {
      ifs_lattice >> neighbor_lattice_id;
      config.third_neighbors_adjacency_list_[lattice_id].push_back(neighbor_lattice_id);
    }
    ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  std::ifstream ifs_element(element_filename, std::ifstream::in);
  if (!ifs_element.is_open()) {
    throw std::runtime_error("Cannot open " + element_filename);
  }
  std::string type;
  config.atom_vector_.reserve(num_atoms);
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id) {
    ifs_element >> type;
    config.atom_vector_.emplace_back(atom_id, type);
    ifs_element.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  std::ifstream ifs_map(map_filename, std::ifstream::in);
  if (!ifs_map.is_open()) {
    throw std::runtime_error("Cannot open " + map_filename);
  }
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
    for (const auto neighbor_lattice_index: first_neighbors_adjacency_list_[i]) {
      ofs << neighbor_lattice_index << ' ';
    }
    for (const auto neighbor_lattice_index: second_neighbors_adjacency_list_[i]) {
      ofs << neighbor_lattice_index << ' ';
    }
    for (const auto neighbor_lattice_index: third_neighbors_adjacency_list_[i]) {
      ofs << neighbor_lattice_index << ' ';
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
  for (size_t atom_id = 0; atom_id < atom_vector_.size(); ++atom_id) {
    ofs << atom_to_lattice_hashmap_.at(atom_id) << std::endl;
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

  const Factor_t num_cells{static_cast<size_t>(std::floor(ScalarLength(basis_[0]) / constants::kNearNeighborsCutoff)),
                           static_cast<size_t>(std::floor(ScalarLength(basis_[1]) / constants::kNearNeighborsCutoff)),
                           static_cast<size_t>(std::floor(ScalarLength(basis_[2]) / constants::kNearNeighborsCutoff))};

  std::vector<std::vector<size_t>> cells{num_cells[0] * num_cells[1] * num_cells[2]};

  for (size_t lattice_id = 0; lattice_id < GetNumAtoms(); ++lattice_id) {
    const auto relative_position = lattice_vector_[lattice_id].GetRelativePosition();
    Factor_t cell_position = ToFactor(ElementFloor(ElementProduct(ToVector(num_cells), relative_position)));
    const auto cell_idx = (cell_position[0] * num_cells[1] + cell_position[1]) * num_cells[2] + cell_position[2];
    cells.at(cell_idx).push_back(lattice_id);
  }

  static const std::vector<std::tuple<int, int, int>> OFFSET_LIST = []() {
    std::vector<std::tuple<int, int, int>> offsets;
    for (int x: {-1, 0, 1}) {
      for (int y: {-1, 0, 1}) {
        for (int z: {-1, 0, 1}) {
          offsets.emplace_back(x, y, z);
        }
      }
    }
    return offsets;
  }();
  // Create neighbor list, iterate over each cell and find neighboring points
  for (size_t cell_idx = 0; cell_idx < cells.size(); ++cell_idx) {
    auto &cell = cells.at(cell_idx);
    const size_t i = cell_idx / (num_cells[1] * num_cells[2]);
    const size_t j = (cell_idx % (num_cells[1] * num_cells[2])) / num_cells[2];
    const size_t k = cell_idx % num_cells[2];
    // Check neighboring cells, taking into account periodic boundaries
    for (auto [di, dj, dk]: OFFSET_LIST) {
      const size_t ni =
          (num_cells[0] + i + (di >= 0 ? static_cast<size_t>(di) : -static_cast<size_t>(-di))) % num_cells[0];
      const size_t nj =
          (num_cells[1] + j + (dj >= 0 ? static_cast<size_t>(dj) : -static_cast<size_t>(-dj))) % num_cells[1];
      const size_t nk =
          (num_cells[2] + k + (dk >= 0 ? static_cast<size_t>(dk) : -static_cast<size_t>(-dk))) % num_cells[2];
      const size_t neighbor_cell_idx = (ni * num_cells[1] + nj) * num_cells[2] + nk;
      auto &neighbor_cell = cells.at(neighbor_cell_idx);
      // For each point in the cell, check if it's close to any point in the neighboring cell
      for (size_t lattice_id1: cell) {
        for (size_t lattice_id2: neighbor_cell) {
          // Make sure we're not comparing a point to itself, and don't double-count pairs within the same cell
          if (lattice_id2 >= lattice_id1) {
            continue;
          }
          // Calculate distance
          Vector_t absolute_distance_vector =
              GetRelativeDistanceVectorLattice(lattice_vector_[lattice_id1], lattice_vector_[lattice_id2]) * basis_;
          if (std::abs(absolute_distance_vector[0]) > constants::kNearNeighborsCutoff) {
            continue;
          }
          if (std::abs(absolute_distance_vector[1]) > constants::kNearNeighborsCutoff) {
            continue;
          }
          if (std::abs(absolute_distance_vector[2]) > constants::kNearNeighborsCutoff) {
            continue;
          }
          const double absolute_distance_square = Inner(absolute_distance_vector);
          // If the distance is less than the cutoff, the points are bonded
          if (absolute_distance_square < first_r_cutoff_square) {
            first_neighbors_adjacency_list_[lattice_id1].push_back(lattice_id2);
            first_neighbors_adjacency_list_[lattice_id2].push_back(lattice_id1);
          } else if (absolute_distance_square < second_r_cutoff_square) {
            second_neighbors_adjacency_list_[lattice_id1].push_back(lattice_id2);
            second_neighbors_adjacency_list_[lattice_id2].push_back(lattice_id1);
          } else if (absolute_distance_square < third_r_cutoff_square) {
            third_neighbors_adjacency_list_[lattice_id1].push_back(lattice_id2);
            third_neighbors_adjacency_list_[lattice_id2].push_back(lattice_id1);
          }
        }
      }
    }
  }

  for (size_t lattice_id = 0; lattice_id < GetNumAtoms(); ++lattice_id) {
    std::sort(first_neighbors_adjacency_list_.at(lattice_id).begin(),
              first_neighbors_adjacency_list_.at(lattice_id).end());
    std::sort(second_neighbors_adjacency_list_.at(lattice_id).begin(),
              second_neighbors_adjacency_list_.at(lattice_id).end());
    std::sort(third_neighbors_adjacency_list_.at(lattice_id).begin(),
              third_neighbors_adjacency_list_.at(lattice_id).end());
  }
}

void RotateLatticeVector(std::vector<Lattice> &lattice_list, const Matrix_t &rotation_matrix) {
  const auto move_distance_after_rotation = Vector_t{0.5, 0.5, 0.5} - (Vector_t{0.5, 0.5, 0.5} * rotation_matrix);
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

Config GenerateFCC(const Factor_t &factors, Element element) {
  Matrix_t basis{{{constants::kLatticeConstant * static_cast<double>(factors[0]), 0, 0},
                  {0, constants::kLatticeConstant * static_cast<double>(factors[1]), 0},
                  {0, 0, constants::kLatticeConstant * static_cast<double>(factors[2])}}};
  const size_t num_atoms = 4 * factors[0] * factors[1] * factors[2];
  auto x_length = static_cast<double>(factors[0]);
  auto y_length = static_cast<double>(factors[1]);
  auto z_length = static_cast<double>(factors[2]);
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(num_atoms);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(num_atoms);
  size_t count = 0;
  for (size_t k = 0; k < factors[2]; ++k) {
    for (size_t j = 0; j < factors[1]; ++j) {
      for (size_t i = 0; i < factors[0]; ++i) {
        auto x_ref = static_cast<double>(i);
        auto y_ref = static_cast<double>(j);
        auto z_ref = static_cast<double>(k);
        std::vector<Vector_t> relative_position_list = {
            {x_ref / x_length, y_ref / y_length, z_ref / z_length},
            {(x_ref + 0.5) / x_length, (y_ref + 0.5) / y_length, z_ref / z_length},
            {(x_ref + 0.5) / x_length, y_ref / y_length, (z_ref + 0.5) / z_length},
            {x_ref / x_length, (y_ref + 0.5) / y_length, (z_ref + 0.5) / z_length}};

        for (const auto &relative_position: relative_position_list) {
          lattice_vector.emplace_back(count, relative_position * basis, relative_position);
          atom_vector.emplace_back(count, element);
          count++;
        }
      }
    }
  }
  return Config{basis, lattice_vector, atom_vector, true};
}

Config GenerateSoluteConfigFromExcitingPure(Config config, const std::map<Element, size_t> &solute_atom_count) {
  std::unordered_set<size_t> unavailable_position{};
  std::unordered_set<size_t> available_position {};
  for (size_t i = 0; i < config.GetNumAtoms(); ++i) {
    available_position.emplace(i);
  }

  static std::mt19937_64 generator(
      static_cast<unsigned long long int>(std::chrono::system_clock::now().time_since_epoch().count()));

  size_t selected_lattice_index{};
  for (const auto &[solute_atom, count]: solute_atom_count) {
    for (size_t it = 0; it < count; ++it) {
      size_t ct = 0;
      do {
        if (ct > 10000) {
          std::cerr << "Size is too small. Cannot generate correct config.\n";
          break;
        }
        ++ct;
        std::vector<size_t> available_positions_vec(available_position.begin(), available_position.end());
        std::uniform_int_distribution<size_t> dis(0, available_positions_vec.size() - 1);
        size_t random_index = dis(generator);
        selected_lattice_index = available_positions_vec[random_index];
      } while (unavailable_position.find(selected_lattice_index) != unavailable_position.end());
      config.SetAtomElementTypeAtLattice(selected_lattice_index, solute_atom);
      available_position.erase(selected_lattice_index);
      unavailable_position.emplace(selected_lattice_index);
      std::copy(config.GetFirstNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config.GetFirstNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position, unavailable_position.begin()));
      std::copy(config.GetSecondNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config.GetSecondNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position, unavailable_position.begin()));
      std::copy(config.GetThirdNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config.GetThirdNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position, unavailable_position.begin()));
    }
  }
  return config;
}

Config GenerateSoluteConfig(const Factor_t &factors,
                            const Element solvent_element,
                            const std::map<Element, size_t> &solute_atom_count) {
  return GenerateSoluteConfigFromExcitingPure(GenerateFCC(factors, solvent_element), solute_atom_count);
}

// Config GenerateClusteredConfigFromExcitingPure(Config config,
//                                                const std::map<Element, size_t> &solute_atom_count) {
//   return Config();
// }
// Config GenerateClusteredConfig(const Factor_t &factors,
//                                Element solvent_element,
//                                const std::map<Element, size_t> &solute_atom_count) {
//   return Config();
// }
}    // namespace cfg
