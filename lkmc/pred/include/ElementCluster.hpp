#ifndef LKMC_LKMC_PRED_INCLUDE_ELEMENTCLUSTER_HPP_
#define LKMC_LKMC_PRED_INCLUDE_ELEMENTCLUSTER_HPP_
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include "Element.hpp"
namespace pred {

class ElementCluster {
  public:
    /// Constructor
    ElementCluster(size_t label, std::vector<Element> element_vector)
        : label_(label), element_vector_(std::move(element_vector)) {
      std::sort(element_vector_.begin(), element_vector_.end());
    }
    template<typename ... Ts>
    explicit ElementCluster(size_t label, Ts &&... ts)
        : label_(label), element_vector_{std::forward<Ts>(ts)...} {
      std::sort(element_vector_.begin(), element_vector_.end());
    }
    /// Getter
    [[nodiscard]] size_t GetSize() const {
      return element_vector_.size();
    }
    [[nodiscard]] size_t GetLabel() const {
      return label_;
    }
    [[nodiscard]] const std::vector<Element> &GetElementVector() const {
      return element_vector_;
    }
    /// Operators
    friend bool operator<(const ElementCluster &lhs, const ElementCluster &rhs) {
      if (lhs.element_vector_.size() < rhs.element_vector_.size()) { return true; }
      if (rhs.element_vector_.size() < lhs.element_vector_.size()) { return false; }
      if (lhs.label_ < rhs.label_) { return true; }
      if (rhs.label_ < lhs.label_) { return false; }
      for (size_t i = 0; i < lhs.element_vector_.size(); ++i) {
        if (lhs.element_vector_[i] < rhs.element_vector_[i]) { return true; }
        if (rhs.element_vector_[i] < lhs.element_vector_[i]) { return false; }
      }
      return false;
    }
    friend bool operator==(const ElementCluster &lhs,
                           const ElementCluster &rhs) {
      if (lhs.element_vector_.size() != rhs.element_vector_.size()) { return false; }
      if (lhs.label_ != rhs.label_) { return false; }
      for (size_t i = 0; i < lhs.element_vector_.size(); ++i) {
        if (lhs.element_vector_[i] != rhs.element_vector_[i])
          return false;
      }
      return true;
    }
    friend std::ostream &operator<<(std::ostream &os, const ElementCluster &element_cluster) {
      os << element_cluster.label_ << ' ';
      for (auto element: element_cluster.element_vector_) {
        os << element.GetString() << '-';
      }
      return os;
    }
    friend size_t hash_value(const ElementCluster &element_cluster) {
      size_t seed = 0;
      boost::hash_combine(seed, element_cluster.element_vector_.size());
      boost::hash_combine(seed, element_cluster.label_);
      for (auto element: element_cluster.element_vector_) {
        boost::hash_combine(seed, element);
      }
      return seed;
    }
  private:
    size_t label_;
    std::vector<Element> element_vector_;
};

inline std::unordered_map<
    ElementCluster, int, boost::hash<ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &type_set) {
  std::unordered_map<ElementCluster, int,
                     boost::hash<ElementCluster> > initialized_cluster_hashmap;

  for (const auto &element1: type_set) {
    initialized_cluster_hashmap[ElementCluster(0, element1)] = 0;
    for (const auto &element2: type_set) {
      if (element2 == Element("X")) {
        continue;
      }
      if (element1 == Element("X") && element2.GetString()[0] == 'p') {
        continue;
      }
      for (size_t label = 1; label <= 3; ++label) {
        initialized_cluster_hashmap[ElementCluster(label, element1, element2)] = 0;
      }
      for (const auto &element3: type_set) {
        if (element3 == Element("X") || element3.GetString()[0] == 'p') {
          continue;
        }
        for (size_t label = 4; label < 11; ++label) {
          initialized_cluster_hashmap[ElementCluster(label, element1, element2, element3)] = 0;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_ELEMENTCLUSTER_HPP_
