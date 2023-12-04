/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 9/23/20 1:29 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 9:50 PM                                                           *
 **************************************************************************************************/

/*! \file  Element.hpp
 *  \brief File for the Element class definition and implementation.
 */

#ifndef LMC_CONSTANT_INCLUDE_ELEMENT_HPP_
#define LMC_CONSTANT_INCLUDE_ELEMENT_HPP_

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>

/*! \enum ElementName
 *  \brief Define the names of the elements that can be used. X is used for the vacancy.
 */
enum class ElementName {
  X,
  H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge,
  As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm,
  Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U,
  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
};

/*! \class Element
 *  \brief A data structure to represent a chemical element.
 *
 *  This class encapsulates the data and operations for a chemical element.
 */
class Element {
 public:
  /*! \brief Default constructor.
   */
  Element() = default;

  /*! \brief Constructor to create an Element instance from an ElementName.
   *  \param element_type : The symbol of this element.
   */
  constexpr explicit Element(ElementName element_type) : element_name_(element_type) {}

  /*! \brief Constructor to create an Element instance from a string.
   *  \param element_string : The name of this element.
   */
  explicit Element(const std::string &element_string) {
    auto it = element_string_map_.find(element_string);
    // If the string does not correspond to a valid element, throw an exception
    if (it == element_string_map_.end()) {
      throw std::invalid_argument("Invalid element string: " + element_string);
    }
    element_name_ = it->second;
  }

  /*! \brief Comparison and other operator overloads.
   */
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
    return GetElementString() < rhs.GetElementString();
  }
  bool operator<(const ElementName rhs) const {
    return GetElementString() < Element(rhs).GetElementString();
  }

  /*! \brief Hash function.
   */
  friend size_t hash_value(Element element) {
    return static_cast<std::size_t>(element.element_name_);
  }

  /*! \brief Returns the atomic index of this element.
   *  \return : The atomic index of this element.
   */
  [[nodiscard]] size_t GetAtomicIndex() const {
    return GetProperties().index;
  }

  /*! \brief Returns the name of this element.
   *  \return : The name of this element.
   */
  [[nodiscard]] std::string GetElementString() const {
    return GetProperties().name;
  }

  /*! \brief Returns the mass of this element.
   *  \return : The mass of this element.
   */
  [[nodiscard]] double GetMass() const {
    return GetProperties().mass;
  }

  /*! \brief stream operator output for this element.
   *  \param os      : The output stream.
   *  \param element : The Element to be streamed.
   *  \return        : The output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const Element &element) {
    os << element.GetElementString();
    return os;
  }

 private:

  /*! \struct ElementProperties
   *  \brief A data structure to hold the properties of an element.
   */
  struct ElementProperties {
    /// The atomic index of this element
    size_t index;
    /// The symbol of this element in string.
    std::string name;
    /// The mass of this element.
    double mass;
  };

  /// Inline static map holding the properties of each ElementName.
  inline static const std::unordered_map<ElementName, ElementProperties> element_properties_map_ = {
      {ElementName::X, {0, "X", 0.0}},
      {ElementName::H, {1, "H", 1.008}}, {ElementName::He, {2, "He", 4.0026}},
      {ElementName::Li, {3, "Li", 6.94}}, {ElementName::Be, {4, "Be", 9.0122}},
      {ElementName::B, {5, "B", 10.81}}, {ElementName::C, {6, "C", 12.011}},
      {ElementName::N, {7, "N", 14.007}}, {ElementName::O, {8, "O", 15.999}},
      {ElementName::F, {9, "F", 18.998}}, {ElementName::Ne, {10, "Ne", 20.180}},
      {ElementName::Na, {11, "Na", 22.990}}, {ElementName::Mg, {12, "Mg", 24.305}},
      {ElementName::Al, {13, "Al", 26.982}}, {ElementName::Si, {14, "Si", 28.085}},
      {ElementName::P, {15, "P", 30.974}}, {ElementName::S, {16, "S", 32.06}},
      {ElementName::Cl, {17, "Cl", 35.45}}, {ElementName::Ar, {18, "Ar", 39.948}},
      {ElementName::K, {19, "K", 39.098}}, {ElementName::Ca, {20, "Ca", 40.078}},
      {ElementName::Sc, {21, "Sc", 44.956}}, {ElementName::Ti, {22, "Ti", 47.867}},
      {ElementName::V, {23, "V", 50.942}}, {ElementName::Cr, {24, "Cr", 51.996}},
      {ElementName::Mn, {25, "Mn", 54.938}}, {ElementName::Fe, {26, "Fe", 55.845}},
      {ElementName::Co, {27, "Co", 58.933}}, {ElementName::Ni, {28, "Ni", 58.693}},
      {ElementName::Cu, {29, "Cu", 63.546}}, {ElementName::Zn, {30, "Zn", 65.38}},
      {ElementName::Ga, {31, "Ga", 69.723}}, {ElementName::Ge, {32, "Ge", 72.630}},
      {ElementName::As, {33, "As", 74.922}}, {ElementName::Se, {34, "Se", 78.971}},
      {ElementName::Br, {35, "Br", 79.904}}, {ElementName::Kr, {36, "Kr", 83.798}},
      {ElementName::Rb, {37, "Rb", 85.468}}, {ElementName::Sr, {38, "Sr", 87.62}},
      {ElementName::Y, {39, "Y", 88.906}}, {ElementName::Zr, {40, "Zr", 91.224}},
      {ElementName::Nb, {41, "Nb", 92.906}}, {ElementName::Mo, {42, "Mo", 95.95}},
      {ElementName::Tc, {43, "Tc", 98}}, {ElementName::Ru, {44, "Ru", 101.07}},
      {ElementName::Rh, {45, "Rh", 102.91}}, {ElementName::Pd, {46, "Pd", 106.42}},
      {ElementName::Ag, {47, "Ag", 107.87}}, {ElementName::Cd, {48, "Cd", 112.41}},
      {ElementName::In, {49, "In", 114.82}}, {ElementName::Sn, {50, "Sn", 118.71}},
      {ElementName::Sb, {51, "Sb", 121.76}}, {ElementName::Te, {52, "Te", 127.60}},
      {ElementName::I, {53, "I", 126.90}}, {ElementName::Xe, {54, "Xe", 131.29}},
      {ElementName::Cs, {55, "Cs", 132.91}}, {ElementName::Ba, {56, "Ba", 137.33}},
      {ElementName::La, {57, "La", 138.91}}, {ElementName::Ce, {58, "Ce", 140.12}},
      {ElementName::Pr, {59, "Pr", 140.91}}, {ElementName::Nd, {60, "Nd", 144.24}},
      {ElementName::Pm, {61, "Pm", 145}}, {ElementName::Sm, {62, "Sm", 150.36}},
      {ElementName::Eu, {63, "Eu", 151.96}}, {ElementName::Gd, {64, "Gd", 157.25}},
      {ElementName::Tb, {65, "Tb", 158.93}}, {ElementName::Dy, {66, "Dy", 162.50}},
      {ElementName::Ho, {67, "Ho", 164.93}}, {ElementName::Er, {68, "Er", 167.26}},
      {ElementName::Tm, {69, "Tm", 168.93}}, {ElementName::Yb, {70, "Yb", 173.05}},
      {ElementName::Lu, {71, "Lu", 174.97}}, {ElementName::Hf, {72, "Hf", 178.49}},
      {ElementName::Ta, {73, "Ta", 180.95}}, {ElementName::W, {74, "W", 183.84}},
      {ElementName::Re, {75, "Re", 186.21}}, {ElementName::Os, {76, "Os", 190.23}},
      {ElementName::Ir, {77, "Ir", 192.22}}, {ElementName::Pt, {78, "Pt", 195.08}},
      {ElementName::Au, {79, "Au", 196.97}}, {ElementName::Hg, {80, "Hg", 200.59}},
      {ElementName::Tl, {81, "Tl", 204.38}}, {ElementName::Pb, {82, "Pb", 207.2}},
      {ElementName::Bi, {83, "Bi", 208.98}}, {ElementName::Po, {84, "Po", 209}},
      {ElementName::At, {85, "At", 210}}, {ElementName::Rn, {86, "Rn", 222}},
      {ElementName::Fr, {87, "Fr", 223}}, {ElementName::Ra, {88, "Ra", 226}},
      {ElementName::Ac, {89, "Ac", 227}}, {ElementName::Th, {90, "Th", 232.04}},
      {ElementName::Pa, {91, "Pa", 231.04}}, {ElementName::U, {92, "U", 238.03}},
      {ElementName::Np, {93, "Np", 237}}, {ElementName::Pu, {94, "Pu", 244}},
      {ElementName::Am, {95, "Am", 243}}, {ElementName::Cm, {96, "Cm", 247}},
      {ElementName::Bk, {97, "Bk", 247}}, {ElementName::Cf, {98, "Cf", 251}},
      {ElementName::Es, {99, "Es", 252}}, {ElementName::Fm, {100, "Fm", 257}},
      {ElementName::Md, {101, "Md", 258}}, {ElementName::No, {102, "No", 259}},
      {ElementName::Lr, {103, "Lr", 266}}, {ElementName::Rf, {104, "Rf", 267}},
      {ElementName::Db, {105, "Db", 270}}, {ElementName::Sg, {106, "Sg", 271}},
      {ElementName::Bh, {107, "Bh", 270}}, {ElementName::Hs, {108, "Hs", 277}},
      {ElementName::Mt, {109, "Mt", 276}}, {ElementName::Ds, {110, "Ds", 281}},
      {ElementName::Rg, {111, "Rg", 282}}, {ElementName::Cn, {112, "Cn", 285}},
      {ElementName::Nh, {113, "Nh", 286}}, {ElementName::Fl, {114, "Fl", 289}},
      {ElementName::Mc, {115, "Mc", 288}}, {ElementName::Lv, {116, "Lv", 293}},
      {ElementName::Ts, {117, "Ts", 294}}, {ElementName::Og, {118, "Og", 294}},
  };

  /// Inline static map to map string names to ElementName.
  inline static const std::unordered_map<std::string, ElementName> element_string_map_ = {
      {"X", ElementName::X},
      {"H", ElementName::H}, {"He", ElementName::He}, {"Li", ElementName::Li}, {"Be", ElementName::Be},
      {"B", ElementName::B}, {"C", ElementName::C}, {"N", ElementName::N}, {"O", ElementName::O},
      {"F", ElementName::F}, {"Ne", ElementName::Ne}, {"Na", ElementName::Na}, {"Mg", ElementName::Mg},
      {"Al", ElementName::Al}, {"Si", ElementName::Si}, {"P", ElementName::P}, {"S", ElementName::S},
      {"Cl", ElementName::Cl}, {"Ar", ElementName::Ar}, {"K", ElementName::K}, {"Ca", ElementName::Ca},
      {"Sc", ElementName::Sc}, {"Ti", ElementName::Ti}, {"V", ElementName::V}, {"Cr", ElementName::Cr},
      {"Mn", ElementName::Mn}, {"Fe", ElementName::Fe}, {"Co", ElementName::Co}, {"Ni", ElementName::Ni},
      {"Cu", ElementName::Cu}, {"Zn", ElementName::Zn}, {"Ga", ElementName::Ga}, {"Ge", ElementName::Ge},
      {"As", ElementName::As}, {"Se", ElementName::Se}, {"Br", ElementName::Br}, {"Kr", ElementName::Kr},
      {"Rb", ElementName::Rb}, {"Sr", ElementName::Sr}, {"Y", ElementName::Y}, {"Zr", ElementName::Zr},
      {"Nb", ElementName::Nb}, {"Mo", ElementName::Mo}, {"Tc", ElementName::Tc}, {"Ru", ElementName::Ru},
      {"Rh", ElementName::Rh}, {"Pd", ElementName::Pd}, {"Ag", ElementName::Ag}, {"Cd", ElementName::Cd},
      {"In", ElementName::In}, {"Sn", ElementName::Sn}, {"Sb", ElementName::Sb}, {"Te", ElementName::Te},
      {"I", ElementName::I}, {"Xe", ElementName::Xe}, {"Cs", ElementName::Cs}, {"Ba", ElementName::Ba},
      {"La", ElementName::La}, {"Ce", ElementName::Ce}, {"Pr", ElementName::Pr}, {"Nd", ElementName::Nd},
      {"Pm", ElementName::Pm}, {"Sm", ElementName::Sm}, {"Eu", ElementName::Eu}, {"Gd", ElementName::Gd},
      {"Tb", ElementName::Tb}, {"Dy", ElementName::Dy}, {"Ho", ElementName::Ho}, {"Er", ElementName::Er},
      {"Tm", ElementName::Tm}, {"Yb", ElementName::Yb}, {"Lu", ElementName::Lu}, {"Hf", ElementName::Hf},
      {"Ta", ElementName::Ta}, {"W", ElementName::W}, {"Re", ElementName::Re}, {"Os", ElementName::Os},
      {"Ir", ElementName::Ir}, {"Pt", ElementName::Pt}, {"Au", ElementName::Au}, {"Hg", ElementName::Hg},
      {"Tl", ElementName::Tl}, {"Pb", ElementName::Pb}, {"Bi", ElementName::Bi}, {"Po", ElementName::Po},
      {"At", ElementName::At}, {"Rn", ElementName::Rn}, {"Fr", ElementName::Fr}, {"Ra", ElementName::Ra},
      {"Ac", ElementName::Ac}, {"Th", ElementName::Th}, {"Pa", ElementName::Pa}, {"U", ElementName::U},
      {"Np", ElementName::Np}, {"Pu", ElementName::Pu}, {"Am", ElementName::Am}, {"Cm", ElementName::Cm},
      {"Bk", ElementName::Bk}, {"Cf", ElementName::Cf}, {"Es", ElementName::Es}, {"Fm", ElementName::Fm},
      {"Md", ElementName::Md}, {"No", ElementName::No}, {"Lr", ElementName::Lr}, {"Rf", ElementName::Rf},
      {"Db", ElementName::Db}, {"Sg", ElementName::Sg}, {"Bh", ElementName::Bh}, {"Hs", ElementName::Hs},
      {"Mt", ElementName::Mt}, {"Ds", ElementName::Ds}, {"Rg", ElementName::Rg}, {"Cn", ElementName::Cn},
      {"Nh", ElementName::Nh}, {"Fl", ElementName::Fl}, {"Mc", ElementName::Mc}, {"Lv", ElementName::Lv},
      {"Ts", ElementName::Ts}, {"Og", ElementName::Og},
  };

  /*! \brief Returns the properties of this element.
   *  \return : The properties of this element.
   */
  [[nodiscard]] const ElementProperties &GetProperties() const {
    auto it = element_properties_map_.find(element_name_);
    if (it == element_properties_map_.end()) {
      throw std::invalid_argument("Unexpected element");
    }
    return it->second;
  }

  /// The symbol of this element.
  ElementName element_name_{};
};

#endif //LMC_CONSTANT_INCLUDE_ELEMENT_HPP_
