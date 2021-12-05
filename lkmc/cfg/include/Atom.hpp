#ifndef LKMC_LKMC_CFG_INCLUDE_ATOM_HPP_
#define LKMC_LKMC_CFG_INCLUDE_ATOM_HPP_
#include <fstream>
#include <map>

#include "Element.hpp"
#include "VectorMatrix.hpp"
namespace cfg {
struct Atom {

  public:
    /// Constructor
    Atom() = default;
    Atom(size_t id, const std::string &element_string)
        : id_(id), element_(element_string) {}
    /// Getter
    [[nodiscard]] size_t GetId() const {
      return id_;
    }
    [[nodiscard]] Element GetElement() const {
      return element_;
    }
    [[nodiscard]] std::string GetElementString() const {
      return element_.GetString();
    }
    [[nodiscard]] double GetMass() const {
      return element_.GetMass();
    }
    /// Setter
  private:
    // atom id
    size_t id_{};
    // element type
    Element element_{};
};
} // namespace cfg

#endif //LKMC_LKMC_CFG_INCLUDE_ATOM_HPP_
