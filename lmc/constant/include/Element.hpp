/**************************************************************************************************
 * Copyright (c) 2021-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 9/23/20 1:29 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 11:56 PM                                                          *
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

 private:

  /*! \struct ElementProperties
   *  \brief A data structure to hold the properties of an element.
   */
  struct ElementProperties {
    /// The symbol of this element in string.
    std::string name;
    /// The mass of this element.
    double mass;
  };

  /// Inline static map holding the properties of each ElementName.
  inline static const std::unordered_map<ElementName, ElementProperties> element_properties_map_ = {
      {ElementName::X, {"X", 0.0}}, {ElementName::H, {"H", 1.008}}, {ElementName::He, {"He", 4.0026}},
      {ElementName::Li, {"Li", 6.94}}, {ElementName::Be, {"Be", 9.0122}}, {ElementName::B, {"B", 10.81}},
      {ElementName::C, {"C", 12.011}}, {ElementName::N, {"N", 14.007}}, {ElementName::O, {"O", 15.999}},
      {ElementName::F, {"F", 18.998}}, {ElementName::Ne, {"Ne", 20.180}}, {ElementName::Na, {"Na", 22.990}},
      {ElementName::Mg, {"Mg", 24.305}}, {ElementName::Al, {"Al", 26.982}}, {ElementName::Si, {"Si", 28.085}},
      {ElementName::P, {"P", 30.974}}, {ElementName::S, {"S", 32.06}}, {ElementName::Cl, {"Cl", 35.45}},
      {ElementName::Ar, {"Ar", 39.948}}, {ElementName::K, {"K", 39.098}}, {ElementName::Ca, {"Ca", 40.078}},
      {ElementName::Sc, {"Sc", 44.956}}, {ElementName::Ti, {"Ti", 47.867}}, {ElementName::V, {"V", 50.942}},
      {ElementName::Cr, {"Cr", 51.996}}, {ElementName::Mn, {"Mn", 54.938}}, {ElementName::Fe, {"Fe", 55.845}},
      {ElementName::Co, {"Co", 58.933}}, {ElementName::Ni, {"Ni", 58.693}}, {ElementName::Cu, {"Cu", 63.546}},
      {ElementName::Zn, {"Zn", 65.38}}, {ElementName::Ga, {"Ga", 69.723}}, {ElementName::Ge, {"Ge", 72.630}},
      {ElementName::As, {"As", 74.922}}, {ElementName::Se, {"Se", 78.971}}, {ElementName::Br, {"Br", 79.904}},
      {ElementName::Kr, {"Kr", 83.798}}, {ElementName::Rb, {"Rb", 85.468}}, {ElementName::Sr, {"Sr", 87.62}},
      {ElementName::Y, {"Y", 88.906}}, {ElementName::Zr, {"Zr", 91.224}}, {ElementName::Nb, {"Nb", 92.906}},
      {ElementName::Mo, {"Mo", 95.95}}, {ElementName::Tc, {"Tc", 98}}, {ElementName::Ru, {"Ru", 101.07}},
      {ElementName::Rh, {"Rh", 102.91}}, {ElementName::Pd, {"Pd", 106.42}}, {ElementName::Ag, {"Ag", 107.87}},
      {ElementName::Cd, {"Cd", 112.41}}, {ElementName::In, {"In", 114.82}}, {ElementName::Sn, {"Sn", 118.71}},
      {ElementName::Sb, {"Sb", 121.76}}, {ElementName::Te, {"Te", 127.60}}, {ElementName::I, {"I", 126.90}},
      {ElementName::Xe, {"Xe", 131.29}}, {ElementName::Cs, {"Cs", 132.91}}, {ElementName::Ba, {"Ba", 137.33}},
      {ElementName::La, {"La", 138.91}}, {ElementName::Ce, {"Ce", 140.12}}, {ElementName::Pr, {"Pr", 140.91}},
      {ElementName::Nd, {"Nd", 144.24}}, {ElementName::Pm, {"Pm", 145}}, {ElementName::Sm, {"Sm", 150.36}},
      {ElementName::Eu, {"Eu", 151.96}}, {ElementName::Gd, {"Gd", 157.25}}, {ElementName::Tb, {"Tb", 158.93}},
      {ElementName::Dy, {"Dy", 162.50}}, {ElementName::Ho, {"Ho", 164.93}}, {ElementName::Er, {"Er", 167.26}},
      {ElementName::Tm, {"Tm", 168.93}}, {ElementName::Yb, {"Yb", 173.05}}, {ElementName::Lu, {"Lu", 174.97}},
      {ElementName::Hf, {"Hf", 178.49}}, {ElementName::Ta, {"Ta", 180.95}}, {ElementName::W, {"W", 183.84}},
      {ElementName::Re, {"Re", 186.21}}, {ElementName::Os, {"Os", 190.23}}, {ElementName::Ir, {"Ir", 192.22}},
      {ElementName::Pt, {"Pt", 195.08}}, {ElementName::Au, {"Au", 196.97}}, {ElementName::Hg, {"Hg", 200.59}},
      {ElementName::Tl, {"Tl", 204.38}}, {ElementName::Pb, {"Pb", 207.2}}, {ElementName::Bi, {"Bi", 208.98}},
      {ElementName::Po, {"Po", 209}}, {ElementName::At, {"At", 210}}, {ElementName::Rn, {"Rn", 222}},
      {ElementName::Fr, {"Fr", 223}}, {ElementName::Ra, {"Ra", 226}}, {ElementName::Ac, {"Ac", 227}},
      {ElementName::Th, {"Th", 232.04}}, {ElementName::Pa, {"Pa", 231.04}}, {ElementName::U, {"U", 238.03}},
      {ElementName::Np, {"Np", 237}}, {ElementName::Pu, {"Pu", 244}}, {ElementName::Am, {"Am", 243}},
      {ElementName::Cm, {"Cm", 247}}, {ElementName::Bk, {"Bk", 247}}, {ElementName::Cf, {"Cf", 251}},
      {ElementName::Es, {"Es", 252}}, {ElementName::Fm, {"Fm", 257}}, {ElementName::Md, {"Md", 258}},
      {ElementName::No, {"No", 259}}, {ElementName::Lr, {"Lr", 266}}, {ElementName::Rf, {"Rf", 267}},
      {ElementName::Db, {"Db", 270}}, {ElementName::Sg, {"Sg", 271}}, {ElementName::Bh, {"Bh", 270}},
      {ElementName::Hs, {"Hs", 277}}, {ElementName::Mt, {"Mt", 276}}, {ElementName::Ds, {"Ds", 281}},
      {ElementName::Rg, {"Rg", 282}}, {ElementName::Cn, {"Cn", 285}}, {ElementName::Nh, {"Nh", 286}},
      {ElementName::Fl, {"Fl", 289}}, {ElementName::Mc, {"Mc", 288}}, {ElementName::Lv, {"Lv", 293}},
      {ElementName::Ts, {"Ts", 294}}, {ElementName::Og, {"Og", 294},
      }};

  /// Inline static map to map string names to ElementName.
  inline static const std::unordered_map<std::string, ElementName> element_string_map_ = {
      {"X", ElementName::X}, {"H", ElementName::H}, {"He", ElementName::He}, {"Li", ElementName::Li},
      {"Be", ElementName::Be}, {"B", ElementName::B}, {"C", ElementName::C}, {"N", ElementName::N},
      {"O", ElementName::O}, {"F", ElementName::F}, {"Ne", ElementName::Ne}, {"Na", ElementName::Na},
      {"Mg", ElementName::Mg}, {"Al", ElementName::Al}, {"Si", ElementName::Si}, {"P", ElementName::P},
      {"S", ElementName::S}, {"Cl", ElementName::Cl}, {"Ar", ElementName::Ar}, {"K", ElementName::K},
      {"Ca", ElementName::Ca}, {"Sc", ElementName::Sc}, {"Ti", ElementName::Ti}, {"V", ElementName::V},
      {"Cr", ElementName::Cr}, {"Mn", ElementName::Mn}, {"Fe", ElementName::Fe}, {"Co", ElementName::Co},
      {"Ni", ElementName::Ni}, {"Cu", ElementName::Cu}, {"Zn", ElementName::Zn}, {"Ga", ElementName::Ga},
      {"Ge", ElementName::Ge}, {"As", ElementName::As}, {"Se", ElementName::Se}, {"Br", ElementName::Br},
      {"Kr", ElementName::Kr}, {"Rb", ElementName::Rb}, {"Sr", ElementName::Sr}, {"Y", ElementName::Y},
      {"Zr", ElementName::Zr}, {"Nb", ElementName::Nb}, {"Mo", ElementName::Mo}, {"Tc", ElementName::Tc},
      {"Ru", ElementName::Ru}, {"Rh", ElementName::Rh}, {"Pd", ElementName::Pd}, {"Ag", ElementName::Ag},
      {"Cd", ElementName::Cd}, {"In", ElementName::In}, {"Sn", ElementName::Sn}, {"Sb", ElementName::Sb},
      {"Te", ElementName::Te}, {"I", ElementName::I}, {"Xe", ElementName::Xe}, {"Cs", ElementName::Cs},
      {"Ba", ElementName::Ba}, {"La", ElementName::La}, {"Ce", ElementName::Ce}, {"Pr", ElementName::Pr},
      {"Nd", ElementName::Nd}, {"Pm", ElementName::Pm}, {"Sm", ElementName::Sm}, {"Eu", ElementName::Eu},
      {"Gd", ElementName::Gd}, {"Tb", ElementName::Tb}, {"Dy", ElementName::Dy}, {"Ho", ElementName::Ho},
      {"Er", ElementName::Er}, {"Tm", ElementName::Tm}, {"Yb", ElementName::Yb}, {"Lu", ElementName::Lu},
      {"Hf", ElementName::Hf}, {"Ta", ElementName::Ta}, {"W", ElementName::W}, {"Re", ElementName::Re},
      {"Os", ElementName::Os}, {"Ir", ElementName::Ir}, {"Pt", ElementName::Pt}, {"Au", ElementName::Au},
      {"Hg", ElementName::Hg}, {"Tl", ElementName::Tl}, {"Pb", ElementName::Pb}, {"Bi", ElementName::Bi},
      {"Po", ElementName::Po}, {"At", ElementName::At}, {"Rn", ElementName::Rn}, {"Fr", ElementName::Fr},
      {"Ra", ElementName::Ra}, {"Ac", ElementName::Ac}, {"Th", ElementName::Th}, {"Pa", ElementName::Pa},
      {"U", ElementName::U}, {"Np", ElementName::Np}, {"Pu", ElementName::Pu}, {"Am", ElementName::Am},
      {"Cm", ElementName::Cm}, {"Bk", ElementName::Bk}, {"Cf", ElementName::Cf}, {"Es", ElementName::Es},
      {"Fm", ElementName::Fm}, {"Md", ElementName::Md}, {"No", ElementName::No}, {"Lr", ElementName::Lr},
      {"Rf", ElementName::Rf}, {"Db", ElementName::Db}, {"Sg", ElementName::Sg}, {"Bh", ElementName::Bh},
      {"Hs", ElementName::Hs}, {"Mt", ElementName::Mt}, {"Ds", ElementName::Ds}, {"Rg", ElementName::Rg},
      {"Cn", ElementName::Cn}, {"Nh", ElementName::Nh}, {"Fl", ElementName::Fl}, {"Mc", ElementName::Mc},
      {"Lv", ElementName::Lv}, {"Ts", ElementName::Ts}, {"Og", ElementName::Og},
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
