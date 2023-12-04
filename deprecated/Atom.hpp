/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/31/20 2:20 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 3:32 PM                                                           *
 **************************************************************************************************/

/*! \file  Atom.hpp
 *  \brief File for the Atom struct definition and implementation.
 */

#ifndef LMC_CONFIG_INCLUDE_ATOM_HPP_
#define LMC_CONFIG_INCLUDE_ATOM_HPP_

#include <fstream>
#include "Element.hpp"

/*! \struct Atom
 *  \brief A data structure to represent an atom.
 *
 *  This structure encapsulates the data and operations for an atom, including
 *  its id, element, and mass.
 */
struct Atom {
 public:
  /*! \brief Default constructor.
   */
  Atom() = default;

  /*! \brief Constructor with id and element object.
   *  \param id      : The unique identifier for this atom.
   *  \param element : The element this atom represents.
   */
  Atom(size_t id, const Element element)
      : id_(id), element_(element) {}

  /*! \brief Constructor with id and element name.
   *  \param id             : The unique identifier for this atom.
   *  \param element_string : The name of the element this atom represents.
   */
  Atom(size_t id, const std::string &element_string)
      : id_(id), element_(element_string) {}

  /*! \brief Returns the id of this atom.
   *  \return : The id of this atom.
   */
  [[nodiscard]] size_t GetId() const {
    return id_;
  }

  /*! \brief Returns the element of this atom.
   *  \return : The element of this atom.
   */
  [[nodiscard]] Element GetElement() const {
    return element_;
  }

  /*! \brief Returns the name of the element this atom represents.
   *  \return : The name of the element.
   */
  [[nodiscard]] std::string GetElementString() const {
    return element_.GetString();
  }

  /*! \brief Returns the mass of the element this atom represents.
   *  \return : The mass of the element.
   */
  [[nodiscard]] double GetMass() const {
    return element_.GetMass();
  }

  /*! \brief Sets the id of this atom.
   *  \param id : The new id of this atom.
   */
  void SetId(size_t id) {
    id_ = id;
  }

  /*! \brief Sets the element of this atom.
   *  \param element : The new element of this atom.
   */
  void SetElement(const Element &element) {
    element_ = element;
  }

 private:
  /// Unique identifier for this atom.
  size_t id_{};

  /// The element this atom represents.
  Element element_{};
};

#endif //LMC_CONFIG_INCLUDE_ATOM_HPP_
