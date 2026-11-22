#include "Element.hpp"

#include <gtest/gtest.h>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct ExpectedElement {
  ElementName name;
  std::string symbol;
  double mass;
};
const std::vector<ExpectedElement> kExpectedElements{
    {ElementName::X, "X", 0.0},       {ElementName::Al, "Al", 26.98},   {ElementName::Mg, "Mg", 24.31},
    {ElementName::Zn, "Zn", 65.38},   {ElementName::Cu, "Cu", 63.55},   {ElementName::Sn, "Sn", 118.71},
    {ElementName::pAl, "pAl", 26.98}, {ElementName::pMg, "pMg", 24.31}, {ElementName::pZn, "pZn", 65.38},
    {ElementName::pCu, "pCu", 63.55}, {ElementName::pSn, "pSn", 118.71}};
}    // namespace

TEST(Element, KnownNamesRoundTrip) {
  for (const std::string name: {"X", "Al", "Mg", "Zn", "Cu", "Sn", "pAl", "pMg", "pZn", "pCu", "pSn"}) {
    EXPECT_EQ(Element(name).GetString(), name);
  }
}

// An unknown name used to become the vacancy "X", which turned a typo into an extra vacancy.
TEST(Element, UnknownNameIsRejected) {
  EXPECT_THROW(Element("Si"), std::invalid_argument);
  EXPECT_THROW(Element("al"), std::invalid_argument);
  EXPECT_THROW(Element(""), std::invalid_argument);
}

// The masses are written to the configuration files and enter the centre of mass of the solutes.
TEST(Element, HasTheNameAndTheMassOfItsSymbol) {
  for (const auto &[name, symbol, mass]: kExpectedElements) {
    const Element element(symbol);
    EXPECT_TRUE(element == name) << symbol;
    EXPECT_TRUE(element == Element(name)) << symbol;
    EXPECT_FALSE(element != name) << symbol;
    EXPECT_EQ(element.GetString(), symbol);
    EXPECT_EQ(element.GetMass(), mass) << symbol;
  }
}

// The number of an element stands for it in the lookup tables of the predictors and in the key of the barrier cache.
TEST(Element, HasTheNumberOfItsName) {
  std::set<size_t> numbers;
  for (const auto &[name, symbol, mass]: kExpectedElements) {
    const Element element(symbol);
    EXPECT_EQ(hash_value(element), static_cast<size_t>(name)) << symbol;
    EXPECT_EQ(static_cast<ElementName>(element), name) << symbol;
    numbers.insert(hash_value(element));
  }
  EXPECT_EQ(numbers.size(), kExpectedElements.size());
}

TEST(Element, DefaultIsTheVacancy) {
  EXPECT_TRUE(Element() == ElementName::X);
}

// The order of the features of the models follows from this order.
TEST(Element, IsOrderedBySymbol) {
  std::set<Element> elements;
  for (const auto &expected: kExpectedElements) {
    elements.emplace(expected.name);
  }
  std::vector<std::string> symbols;
  for (const auto element: elements) {
    symbols.push_back(element.GetString());
  }
  EXPECT_EQ(symbols,
            (std::vector<std::string>{"Al", "Cu", "Mg", "Sn", "X", "Zn", "pAl", "pCu", "pMg", "pSn", "pZn"}));
  EXPECT_TRUE(Element("Al") < ElementName::Cu);
  EXPECT_FALSE(Element("Cu") < ElementName::Al);
  EXPECT_FALSE(Element("Al") < ElementName::Al);
}

TEST(Element, PseudoElementHasTheMassOfItsElement) {
  for (const std::string symbol: {"Al", "Mg", "Zn", "Cu", "Sn"}) {
    const auto pseudo = Element(symbol).GetPseudo();
    EXPECT_EQ(pseudo.GetString(), "p" + symbol);
    EXPECT_EQ(pseudo.GetMass(), Element(symbol).GetMass()) << symbol;
  }
}

TEST(Element, VacancyAndPseudoElementsHaveNoPseudoElement) {
  EXPECT_THROW(static_cast<void>(Element("X").GetPseudo()), std::invalid_argument);
  for (const std::string symbol: {"pAl", "pMg", "pZn", "pCu", "pSn"}) {
    EXPECT_THROW(static_cast<void>(Element(symbol).GetPseudo()), std::invalid_argument) << symbol;
  }
}

TEST(Element, NumberOutsideTheEnumerationIsRejected) {
  const Element element(static_cast<ElementName>(1000));
  EXPECT_THROW(static_cast<void>(element.GetString()), std::invalid_argument);
  EXPECT_THROW(static_cast<void>(element.GetMass()), std::invalid_argument);
  EXPECT_THROW(static_cast<void>(element.GetPseudo()), std::invalid_argument);
}
