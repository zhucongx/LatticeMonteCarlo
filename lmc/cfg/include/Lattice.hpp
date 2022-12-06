#ifndef LMC_LMC_CFG_INCLUDE_LATTICE_HPP_
#define LMC_LMC_CFG_INCLUDE_LATTICE_HPP_
#include "Constants.hpp"
#include "VectorMatrix.hpp"
namespace cfg {

class Lattice {
  public:
    /// Constructor
    Lattice() = default;
    Lattice(size_t id, const Vector_t &position)
        : id_(id),
          cartesian_position_(position),
          relative_position_(position) {}
    Lattice(size_t id,
            const Vector_t &cartesian_position,
            const Vector_t &relative_position)
        : id_(id),
          cartesian_position_(cartesian_position),
          relative_position_(relative_position) {}
    Lattice(size_t id, double x, double y, double z)
        : id_(id),
          cartesian_position_{x, y, z},
          relative_position_{x, y, z} {}
    /// Getter
    [[nodiscard]] size_t GetId() const {
      return id_;
    }
    [[nodiscard]] const Vector_t &GetCartesianPosition() const {
      return cartesian_position_;
    }
    [[nodiscard]] const Vector_t &GetRelativePosition() const {
      return relative_position_;
    }
    /// Setter
    void SetId(size_t id) {
      id_ = id;
    }
    void SetCartesianPosition(const Vector_t &cartesian_position) {
      cartesian_position_ = cartesian_position;
    }
    void SetRelativePosition(const Vector_t &relative_position) {
      relative_position_ = relative_position;
    }
  private:
    // lattice id
    size_t id_{};
    // absolute position
    Vector_t cartesian_position_{};
    // relative position in the box
    Vector_t relative_position_{};
};

inline Vector_t GetRelativeDistanceVectorLattice(const Lattice &first, const Lattice &second) {
  Vector_t relative_distance_vector = second.GetRelativePosition() - first.GetRelativePosition();
  // periodic boundary conditions
  for (const auto kDim: All_Dimensions) {
    while (relative_distance_vector[kDim] >= 0.5)
      relative_distance_vector[kDim] -= 1;
    while (relative_distance_vector[kDim] < -0.5)
      relative_distance_vector[kDim] += 1;
  }
  return relative_distance_vector;
}
} // cfg
#endif //LMC_LMC_CFG_INCLUDE_LATTICE_HPP_
