#ifndef LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
#define LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>

// The elements, one per line: symbol and mass. A new element is one new line. Its name, its mass, its lookup by symbol
// and its pseudo element, which has the symbol with a "p" in front and the same mass, follow from this list. The
// vacancy X is not in the list: it has no pseudo element.
#define LMC_ELEMENT_LIST(ELEMENT) \
  ELEMENT(Al, 26.98)              \
  ELEMENT(Mg, 24.31)              \
  ELEMENT(Zn, 65.38)              \
  ELEMENT(Cu, 63.55)              \
  ELEMENT(Sn, 118.71)             \
  /* A new element goes above this line. */

#define LMC_ELEMENT_NAME(symbol, mass) symbol,
#define LMC_PSEUDO_ELEMENT_NAME(symbol, mass) p##symbol,
// The vacancy, the elements of the list and their pseudo elements.
enum class ElementName { X, LMC_ELEMENT_LIST(LMC_ELEMENT_NAME) LMC_ELEMENT_LIST(LMC_PSEUDO_ELEMENT_NAME) };
#undef LMC_ELEMENT_NAME
#undef LMC_PSEUDO_ELEMENT_NAME

class Element {
 public:
  Element() = default;
  constexpr explicit Element(ElementName element_type)
      : element_name_(element_type) {}
  explicit Element(const std::string &element_string) {
    for (size_t number = 0; number < std::size(kElementData); ++number) {
      if (kElementData[number].symbol == element_string) {
        element_name_ = static_cast<ElementName>(number);
        return;
      }
    }
    // Mapping an unknown name to the vacancy "X" would silently turn a typo into extra vacancies.
    throw std::invalid_argument("Unknown element: " + element_string);
  }
  constexpr explicit operator ElementName() const {
    return element_name_;
  }
  explicit operator bool() = delete;
  constexpr bool operator==(Element rhs) const {
    return element_name_ == rhs.element_name_;
  }
  constexpr bool operator!=(Element rhs) const {
    return element_name_ != rhs.element_name_;
  }
  constexpr bool operator==(ElementName rhs) const {
    return element_name_ == rhs;
  }
  constexpr bool operator!=(ElementName rhs) const {
    return element_name_ != rhs;
  }
  // By symbol. The order of the features of the models follows from this order.
  bool operator<(Element rhs) const {
    return GetData().symbol < rhs.GetData().symbol;
  }
  bool operator<(const ElementName rhs) const {
    return *this < Element(rhs);
  }
  friend size_t hash_value(Element element) {
    return static_cast<std::size_t>(element.element_name_);
  }
  [[nodiscard]] std::string GetString() const {
    return std::string(GetData().symbol);
  }
  [[nodiscard]] double GetMass() const {
    return GetData().mass;
  }
  [[nodiscard]] Element GetPseudo() const {
    // The vacancy has number 0, and the numbers of the pseudo elements follow those of the elements.
    const auto number = static_cast<size_t>(element_name_);
    if (number == 0 || number > kNumElements) {
      throw std::invalid_argument("Unexpected pseudo element");
    }
    return Element(static_cast<ElementName>(number + kNumElements));
  }

 private:
  struct Data {
    std::string_view symbol;
    double mass;
  };
#define LMC_ELEMENT_DATA(symbol, mass) {#symbol, mass},
#define LMC_PSEUDO_ELEMENT_DATA(symbol, mass) {"p" #symbol, mass},
  // In the order of ElementName.
  static constexpr Data kElementData[]{
      {"X", 0.00}, LMC_ELEMENT_LIST(LMC_ELEMENT_DATA) LMC_ELEMENT_LIST(LMC_PSEUDO_ELEMENT_DATA)};
#undef LMC_ELEMENT_DATA
#undef LMC_PSEUDO_ELEMENT_DATA
  // Without the vacancy and the pseudo elements.
  static constexpr size_t kNumElements = (std::size(kElementData) - 1) / 2;
  // The predictors use the number of an element as one byte.
  static_assert(std::size(kElementData) <= 256, "The number of an element has to fit into one byte");

  [[nodiscard]] const Data &GetData() const {
    const auto number = static_cast<size_t>(element_name_);
    if (number >= std::size(kElementData)) {
      throw std::invalid_argument("Unexpected element");
    }
    return kElementData[number];
  }

  ElementName element_name_{};
};
#undef LMC_ELEMENT_LIST
#endif    //LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
