#ifndef LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
#define LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
#include <cmath>
#include <string>
#include <map>
#include <iostream>

enum class ElementName {
  X, Al, Mg, Zn, Cu, Sn, pAl, pMg, pZn, pCu, pSn
};

class Element {
  public:
    Element() = default;
    constexpr explicit Element(ElementName element_type) : element_name_(element_type) {}
    explicit Element(const std::string &element_string) {
      const std::map<std::string, ElementName> ElementStrings{
          {"X", ElementName::X},
          {"Al", ElementName::Al},
          {"Mg", ElementName::Mg},
          {"Zn", ElementName::Zn},
          {"Cu", ElementName::Cu},
          {"Sn", ElementName::Sn},
          {"pAl", ElementName::pAl},
          {"pMg", ElementName::pMg},
          {"pZn", ElementName::pZn},
          {"pCu", ElementName::pCu},
          {"pSn", ElementName::pSn}
      };
      auto it = ElementStrings.find(element_string);
      element_name_ = it == ElementStrings.end() ? ElementName::X : it->second;
    }

    constexpr explicit operator ElementName() const { return element_name_; }
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
    bool operator<(Element rhs) const {
      return GetString() < rhs.GetString();
    }
    bool operator<(const ElementName rhs) const {
      return GetString() < Element(rhs).GetString();
    }
    friend size_t hash_value(Element element) {
      return static_cast<std::size_t>(element.element_name_);
    }
    [[nodiscard]] std::string GetString() const {
      switch (element_name_) {
        case ElementName::X : return "X";
        case ElementName::Al: return "Al";
        case ElementName::Mg: return "Mg";
        case ElementName::Zn: return "Zn";
        case ElementName::Cu: return "Cu";
        case ElementName::Sn: return "Sn";
        case ElementName::pAl: return "pAl";
        case ElementName::pMg: return "pMg";
        case ElementName::pZn: return "pZn";
        case ElementName::pCu: return "pCu";
        case ElementName::pSn: return "pSn";
        default: throw std::invalid_argument("Unexpected pseudo element");
          // omit default case to trigger compiler warning for missing cases
      }
    }
    [[nodiscard]] double GetMass() const {
      switch (element_name_) {
        case ElementName::X : return 0.00;
        case ElementName::Al: return 26.98;
        case ElementName::Mg: return 24.31;
        case ElementName::Zn: return 65.38;
        case ElementName::Cu: return 63.55;
        case ElementName::Sn: return 118.71;
        case ElementName::pAl: return 26.98;
        case ElementName::pMg: return 24.31;
        case ElementName::pZn: return 65.38;
        case ElementName::pCu: return 63.55;
        case ElementName::pSn: return 118.71;
        default: throw std::invalid_argument("Unexpected pseudo element");
          // omit default case to trigger compiler warning for missing cases
      }
    }
    [[nodiscard]] Element GetPseudo() const {
      switch (element_name_) {
        case ElementName::Al : return Element(ElementName::pAl);
        case ElementName::Mg: return Element(ElementName::pMg);
        case ElementName::Zn: return Element(ElementName::pZn);
        case ElementName::Cu: return Element(ElementName::pCu);
        case ElementName::Sn: return Element(ElementName::pSn);
          // omit default case to trigger compiler warning for missing cases
        default: throw std::invalid_argument("Unexpected pseudo element");
      }
    }
  private:
    ElementName element_name_{};
};

#endif //LMC_LMC_CFG_INCLUDE_ELEMENT_HPP_
