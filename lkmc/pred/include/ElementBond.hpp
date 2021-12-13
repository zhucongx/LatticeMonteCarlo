#ifndef LKMC_LKMC_PRED_INCLUDE_ELEMENTBOND_HPP_
#define LKMC_LKMC_PRED_INCLUDE_ELEMENTBOND_HPP_
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
class ElementBond {
  public:
    /// Constructor
    ElementBond(size_t label, const Element &type1, const Element &type2)
        : label_(label), type1_(type1), type2_(type2) {
      if (type2_ < type1_)
        std::swap(type1_, type2_);
    }
    /// Getter
    [[nodiscard]] const Element &GetType1() const {
      return type1_;
    }
    [[nodiscard]] const Element &GetType2() const {
      return type2_;
    }
    /// Operators
    friend bool operator<(const ElementBond &lhs, const ElementBond &rhs) {
      if (lhs.label_ < rhs.label_)
        return true;
      if (rhs.label_ < lhs.label_)
        return false;
      if (lhs.type1_ < rhs.type1_)
        return true;
      if (rhs.type1_ < lhs.type1_)
        return false;
      return lhs.type2_ < rhs.type2_;
    }
    friend bool operator==(const ElementBond &lhs, const ElementBond &rhs) {
      return lhs.label_ == rhs.label_ && lhs.type1_ == rhs.type1_ && lhs.type2_ == rhs.type2_;
    }
    friend size_t hash_value(const ElementBond &bond) {
      size_t seed = 0;
      boost::hash_combine(seed, bond.label_);
      boost::hash_combine(seed, bond.type1_);
      boost::hash_combine(seed, bond.type2_);
      return seed;
    }
  private:
    size_t label_;
    Element type1_;
    Element type2_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ELEMENTBOND_HPP_
