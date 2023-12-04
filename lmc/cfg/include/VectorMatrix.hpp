#ifndef LMC_LMC_CFG_INCLUDE_VECTORMATRIX_HPP_
#define LMC_LMC_CFG_INCLUDE_VECTORMATRIX_HPP_

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "Constants.hpp"
constexpr size_t kDimension = 3;

enum Dimension { kXDimension, kYDimension, kZDimension };
constexpr std::array<Dimension, 3> All_Dimensions{kXDimension, kYDimension, kZDimension};

// By default, it is always a 1 by 3 vector
using Vector_t = std::array<double, kDimension>;
using Matrix_t = std::array<Vector_t, kDimension>;
using Factor_t = std::array<size_t, kDimension>;

inline Factor_t ToFactor(const Vector_t &vector)
{
  return {static_cast<size_t>(vector[0]), static_cast<size_t>(vector[1]), static_cast<size_t>(vector[2])};
}

inline Vector_t ToVector(const Factor_t &factor)
{
  return {static_cast<double>(factor[0]), static_cast<double>(factor[1]), static_cast<double>(factor[2])};
}

inline std::ostream &operator<<(std::ostream &os, const Vector_t &vector)
{
  os << std::fixed << vector[0] << ' ' << vector[1] << ' ' << vector[2];
  return os;
}

inline std::istream &operator>>(std::istream &is, Vector_t &vector)
{
  is >> vector[0] >> vector[1] >> vector[2];
  return is;
}

constexpr double kEpsilon = 1e-8;

inline bool operator==(const Vector_t &lhs, const Vector_t &rhs)
{
  return std::abs(lhs[0] - rhs[0]) < kEpsilon && std::abs(lhs[1] - rhs[1]) < kEpsilon &&
      std::abs(lhs[2] - rhs[2]) < kEpsilon;
}

inline bool operator!=(const Vector_t &lhs, const Vector_t &rhs)
{
  return !(rhs == lhs);
}

inline bool operator<(const Vector_t &lhs, const Vector_t &rhs)
{
  const double x_diff = lhs[0] - rhs[0];
  if (x_diff < -kEpsilon) return true;
  if (x_diff > kEpsilon) return false;
  const double y_diff = lhs[1] - rhs[1];
  if (y_diff < -kEpsilon) return true;
  if (y_diff > kEpsilon) return false;

  return lhs[2] < rhs[2] - kEpsilon;
}

inline Vector_t operator-(const Vector_t &vector)
{
  return {-vector[0], -vector[1], -vector[2]};
}

inline Vector_t &operator+=(Vector_t &lhs, const Vector_t &rhs)
{
  lhs[0] += rhs[0];
  lhs[1] += rhs[1];
  lhs[2] += rhs[2];
  return lhs;
}

inline Vector_t &operator-=(Vector_t &lhs, const Vector_t &rhs)
{
  lhs[0] -= rhs[0];
  lhs[1] -= rhs[1];
  lhs[2] -= rhs[2];
  return lhs;
}

inline Vector_t &operator*=(Vector_t &lhs, double factor)
{
  lhs[0] *= factor;
  lhs[1] *= factor;
  lhs[2] *= factor;
  return lhs;
}

inline Vector_t &operator/=(Vector_t &lhs, double divisor)
{
  lhs[0] /= divisor;
  lhs[1] /= divisor;
  lhs[2] /= divisor;
  return lhs;
}

inline Vector_t operator+(const Vector_t &lhs, const Vector_t &rhs)
{
  Vector_t temp(lhs);
  return (temp += rhs);
}

inline Vector_t operator-(const Vector_t &lhs, const Vector_t &rhs)
{
  Vector_t temp(lhs);
  return (temp -= rhs);
}

inline Vector_t operator+(const Vector_t &lhs, double rhs)
{
  Vector_t temp{rhs, rhs, rhs};
  return (lhs + temp);
}

inline Vector_t operator-(const Vector_t &lhs, double rhs)
{
  Vector_t temp{rhs, rhs, rhs};
  return (lhs - temp);
}

inline Vector_t operator*(const Vector_t &vector, double factor)
{
  Vector_t temp(vector);
  return (temp *= factor);
}

inline Vector_t operator*(double factor, const Vector_t &vector)
{
  return operator*(vector, factor);
}

inline Vector_t operator/(const Vector_t &vector, double divisor)
{
  Vector_t temp(vector);
  return (temp /= divisor);
}

inline double Max(const Vector_t &vector)
{
  return std::max(std::max(vector[0], vector[1]), vector[2]);
}

inline double Min(const Vector_t &vector)
{
  return std::min(std::min(vector[0], vector[1]), vector[2]);
}

inline double Sum(const Vector_t &vector)
{
  return vector[0] + vector[1] + vector[2];
}

inline Vector_t ElementAbs(const Vector_t &vector)
{
  return {std::abs(vector[0]), std::abs(vector[1]), std::abs(vector[2])};
}

inline Vector_t ElementFloor(const Vector_t &vector)
{
  return {std::floor(vector[0]), std::floor(vector[1]), std::floor(vector[2])};
}

inline Vector_t Cross(const Vector_t &first, const Vector_t &second)
{
  return {first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2],
          first[0] * second[1] - first[1] * second[0]};
}

inline double Dot(const Vector_t &first, const Vector_t &second)
{
  return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

inline double Inner(const Vector_t &vector)
{
  return vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
}

inline Vector_t ElementProduct(const Vector_t &first, const Vector_t &second)
{
  return {first[0] * second[0], first[1] * second[1], first[2] * second[2]};
}

inline Vector_t ElementDivide(const Vector_t &dividend, const Vector_t &divisor)
{
  return {dividend[0] / divisor[0], dividend[1] / divisor[1], dividend[2] / divisor[2]};
}

inline double ScalarLength(const Vector_t &vector)
{
  return std::sqrt(Inner(vector));
}

inline Vector_t Normalize(const Vector_t &vector)
{
  double factor = 1.0 / ScalarLength(vector);
  return vector * factor;
}

inline std::ostream &operator<<(std::ostream &os, const Matrix_t &matrix)
{
  os << matrix[0] << '\n' << matrix[1] << '\n' << matrix[2];
  return os;
}

inline std::istream &operator>>(std::istream &is, Matrix_t &matrix)
{
  is >> matrix[0] >> matrix[1] >> matrix[2];
  return is;
}

inline Matrix_t &operator*=(Matrix_t &lhs, double factor)
{
  lhs[0] *= factor;
  lhs[1] *= factor;
  lhs[2] *= factor;
  return lhs;
}

inline Matrix_t &operator/=(Matrix_t &lhs, double divisor)
{
  lhs[0] /= divisor;
  lhs[1] /= divisor;
  lhs[2] /= divisor;
  return lhs;
}

inline Matrix_t operator*(const Matrix_t &matrix, double factor)
{
  Matrix_t temp(matrix);
  return (temp *= factor);
}

inline Matrix_t operator*(double factor, const Matrix_t &matrix)
{
  return operator*(matrix, factor);
}

inline Matrix_t operator/(const Matrix_t &matrix, double divisor)
{
  Matrix_t temp(matrix);
  return (temp /= divisor);
}

inline Vector_t operator*(const Vector_t &lhs, const Matrix_t &rhs)
{
  return {lhs[0] * rhs[0][0] + lhs[1] * rhs[1][0] + lhs[2] * rhs[2][0],
          lhs[0] * rhs[0][1] + lhs[1] * rhs[1][1] + lhs[2] * rhs[2][1],
          lhs[0] * rhs[0][2] + lhs[1] * rhs[1][2] + lhs[2] * rhs[2][2]};
}

inline double Determinant(const Matrix_t &input)
{
  return (input[0][0] * input[1][1] * input[2][2] - input[0][0] * input[1][2] * input[2][1] -
          input[0][1] * input[1][0] * input[2][2] + input[0][1] * input[1][2] * input[2][0] +
          input[0][2] * input[1][0] * input[2][1] - input[0][2] * input[1][1] * input[2][0]);
}

inline Matrix_t TransposeMatrix(const Matrix_t &input)
{
  return {{{input[0][0], input[1][0], input[2][0]},
           {input[0][1], input[1][1], input[2][1]},
           {input[0][2], input[1][2], input[2][2]}}};
}

inline Matrix_t InverseMatrix(const Matrix_t &input)
{
  double det = Determinant(input);
  return {{{(input[1][1] * input[2][2] - input[1][2] * input[2][1]) / det,
            (input[0][2] * input[2][1] - input[0][1] * input[2][2]) / det,
            (input[0][1] * input[1][2] - input[0][2] * input[1][1]) / det},
           {(input[1][2] * input[2][0] - input[1][0] * input[2][2]) / det,
            (input[0][0] * input[2][2] - input[0][2] * input[2][0]) / det,
            (input[0][2] * input[1][0] - input[0][0] * input[1][2]) / det},
           {(input[1][0] * input[2][1] - input[1][1] * input[2][0]) / det,
            (input[0][1] * input[2][0] - input[0][0] * input[2][1]) / det,
            (input[0][0] * input[1][1] - input[0][1] * input[1][0]) / det}}};
}
#endif    //LMC_LMC_CFG_INCLUDE_VECTORMATRIX_HPP_
