#ifndef LKMC_LKMC_CFG_INCLUDE_ELEMENT_HPP_
#define LKMC_LKMC_CFG_INCLUDE_ELEMENT_HPP_
#include <cmath>
#include <string>
#include <map>
#include <iostream>

enum class ElementType {
  X, Al, Mg, Zn, Cu, Sn, pAl, pMg, pZn, pCu, pSn
};
class Element {
  public:
    Element() = default;
    constexpr explicit Element(ElementType element_type) : element_type_(element_type) {}
    explicit Element(const std::string &element_string) {
      const std::map<std::string, ElementType> ElementStrings{
          {"X", ElementType::X},
          {"Al", ElementType::Al},
          {"Mg", ElementType::Mg},
          {"Zn", ElementType::Zn},
          {"Cu", ElementType::Cu},
          {"Sn", ElementType::Sn},
          {"pAl", ElementType::pAl},
          {"pMg", ElementType::pMg},
          {"pZn", ElementType::pZn},
          {"pCu", ElementType::pCu},
          {"pSn", ElementType::pSn}
      };
      auto it = ElementStrings.find(element_string);
      element_type_ = it == ElementStrings.end() ? ElementType::X : it->second;
    }

    constexpr explicit operator ElementType() const { return element_type_; }
    explicit operator bool() = delete;
    constexpr bool operator==(Element rhs) const {
      return element_type_ == rhs.element_type_;
    }
    constexpr bool operator!=(Element rhs) const {
      return element_type_ != rhs.element_type_;
    }
    constexpr bool operator==(ElementType rhs) const {
      return element_type_ == rhs;
    }
    constexpr bool operator!=(ElementType rhs) const {
      return element_type_ != rhs;
    }
    bool operator<(Element rhs) const {
      return GetString() < rhs.GetString();
    }
    bool operator<(const ElementType rhs) const {
      return GetString() < Element(rhs).GetString();
    }
    friend size_t hash_value(Element element) {
      return static_cast<std::size_t>(element.element_type_);
    }
    [[nodiscard]] std::string GetString() const {
      switch (element_type_) {
        case ElementType::X : return "X";
        case ElementType::Al: return "Al";
        case ElementType::Mg: return "Mg";
        case ElementType::Zn: return "Zn";
        case ElementType::Cu: return "Cu";
        case ElementType::Sn: return "Sn";
        case ElementType::pAl: return "pAl";
        case ElementType::pMg: return "pMg";
        case ElementType::pZn: return "pZn";
        case ElementType::pCu: return "pCu";
        case ElementType::pSn: return "pSn";
        default: std::cerr << "Unexpected pseudo element" << std::endl;
          return "undefined element";
          // omit default case to trigger compiler warning for missing cases
      }
    }
    [[nodiscard]] double GetMass() const {
      switch (element_type_) {
        case ElementType::X : return 0.00;
        case ElementType::Al: return 26.98;
        case ElementType::Mg: return 24.31;
        case ElementType::Zn: return 65.38;
        case ElementType::Cu: return 63.55;
        case ElementType::Sn: return 118.71;
        case ElementType::pAl: return 26.98;
        case ElementType::pMg: return 24.31;
        case ElementType::pZn: return 65.38;
        case ElementType::pCu: return 63.55;
        case ElementType::pSn: return 118.71;
        default: std::cerr << "Unexpected pseudo element" << std::endl;
          return 0;
          // omit default case to trigger compiler warning for missing cases
      }
    }
    [[nodiscard]] Element GetPseudo() const {
      switch (element_type_) {
        case ElementType::Al : return Element(ElementType::pAl);
        case ElementType::Mg: return Element(ElementType::pMg);
        case ElementType::Zn: return Element(ElementType::pZn);
        case ElementType::Cu: return Element(ElementType::pCu);
        case ElementType::Sn: return Element(ElementType::pSn);
          // omit default case to trigger compiler warning for missing cases
        default: std::cerr << "Unexpected pseudo element" << std::endl;
          return Element(ElementType::X);
      }
    }
  private:
    ElementType element_type_;
};

#endif //LKMC_LKMC_CFG_INCLUDE_ELEMENT_HPP_