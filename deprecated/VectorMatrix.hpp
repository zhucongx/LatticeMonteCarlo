/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 05/22/20 2:12 AM                                                                        *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/2/23 11:14 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_CFG_INCLUDE_VECTORMATRIX_HPP_
#define LMC_CFG_INCLUDE_VECTORMATRIX_HPP_

#include <cmath>
#include <array>
#include <iostream>
#include <numeric>
#include <iomanip>

#include "Constants.hpp"
constexpr size_t kDimension = 3;

enum Dimension { kXDimension, kYDimension, kZDimension };
constexpr std::array<Dimension, 3>
    All_Dimensions{kXDimension, kYDimension, kZDimension};


// By default, it is always a 1 by 3 vector
using Vector_t = std::array<double, kDimension>;
using Matrix_t = std::array<Vector_t, kDimension>;
using Factor_t = std::array<size_t, kDimension>;

inline std::ostream &operator<<(std::ostream &os, const Vector_t &vector) {
  os << std::fixed << vector[kXDimension] << ' ' << vector[kYDimension] << ' '
     << vector[kZDimension];
  return os;
}

inline std::istream &operator>>(std::istream &is, Vector_t &vector) {
  is >> vector[kXDimension] >> vector[kYDimension] >> vector[kZDimension];
  return is;
}

constexpr double kEpsilon = 1e-8;

inline bool operator==(const Vector_t &lhs, const Vector_t &rhs) {
  return std::abs(lhs[kXDimension] - rhs[kXDimension]) < kEpsilon &&
      std::abs(lhs[kYDimension] - rhs[kYDimension]) < kEpsilon &&
      std::abs(lhs[kZDimension] - rhs[kZDimension]) < kEpsilon;
}

inline bool operator!=(const Vector_t &lhs, const Vector_t &rhs) {
  return !(rhs == lhs);
}

inline bool operator<(const Vector_t &lhs, const Vector_t &rhs) {
  const double x_diff = lhs[kXDimension] - rhs[kXDimension];
  if (x_diff < -kEpsilon)
    return true;
  if (x_diff > kEpsilon)
    return false;
  const double y_diff = lhs[kYDimension] - rhs[kYDimension];
  if (y_diff < -kEpsilon)
    return true;
  if (y_diff > kEpsilon)
    return false;

  return lhs[kZDimension] < rhs[kZDimension] - kEpsilon;
}

inline Vector_t operator-(const Vector_t &vector) {
  return {-vector[kXDimension],
          -vector[kYDimension],
          -vector[kZDimension]};
}

inline Vector_t &operator+=(Vector_t &lhs, const Vector_t &rhs) {
  lhs[kXDimension] += rhs[kXDimension];
  lhs[kYDimension] += rhs[kYDimension];
  lhs[kZDimension] += rhs[kZDimension];
  return lhs;
}

inline Vector_t &operator-=(Vector_t &lhs, const Vector_t &rhs) {
  lhs[kXDimension] -= rhs[kXDimension];
  lhs[kYDimension] -= rhs[kYDimension];
  lhs[kZDimension] -= rhs[kZDimension];
  return lhs;
}

inline Vector_t &operator*=(Vector_t &lhs, double factor) {
  lhs[kXDimension] *= factor;
  lhs[kYDimension] *= factor;
  lhs[kZDimension] *= factor;
  return lhs;
}

inline Vector_t &operator/=(Vector_t &lhs, double divisor) {
  lhs[kXDimension] /= divisor;
  lhs[kYDimension] /= divisor;
  lhs[kZDimension] /= divisor;
  return lhs;
}

inline Vector_t operator+(const Vector_t &lhs, const Vector_t &rhs) {
  Vector_t temp(lhs);
  return (temp += rhs);
}

inline Vector_t operator-(const Vector_t &lhs, const Vector_t &rhs) {
  Vector_t temp(lhs);
  return (temp -= rhs);
}
inline Vector_t operator+(const Vector_t &lhs, double rhs) {
  Vector_t temp{rhs, rhs, rhs};
  return (lhs + temp);
}
inline Vector_t operator-(const Vector_t &lhs, double rhs) {
  Vector_t temp{rhs, rhs, rhs};
  return (lhs - temp);
}

inline Vector_t operator*(const Vector_t &vector, double factor) {
  Vector_t temp(vector);
  return (temp *= factor);
}

inline Vector_t operator*(double factor, const Vector_t &vector) {
  return operator*(vector, factor);
}

inline Vector_t operator/(const Vector_t &vector, double divisor) {
  Vector_t temp(vector);
  return (temp /= divisor);
}

inline double Max(const Vector_t &vector) {
  return std::max(std::max(vector[kXDimension], vector[kYDimension]),
                  vector[kZDimension]);
}

inline double Min(const Vector_t &vector) {
  return std::min(std::min(vector[kXDimension], vector[kYDimension]),
                  vector[kZDimension]);
}

inline double Sum(const Vector_t &vector) {
  return vector[kXDimension] + vector[kYDimension] + vector[kZDimension];
}

inline Vector_t ElementAbs(const Vector_t &vector) {
  return {
      std::abs(vector[kXDimension]), std::abs(vector[kYDimension]),
      std::abs(vector[kZDimension])
  };
}

inline Vector_t ElementFloor(const Vector_t &vector) {
  return {
      std::floor(vector[kXDimension]), std::floor(vector[kYDimension]),
      std::floor(vector[kZDimension])
  };
}

inline Vector_t Cross(const Vector_t &first, const Vector_t &second) {
  return {
      first[kYDimension] * second[kZDimension]
          - first[kZDimension] * second[kYDimension],
      first[kZDimension] * second[kXDimension]
          - first[kXDimension] * second[kZDimension],
      first[kXDimension] * second[kYDimension]
          - first[kYDimension] * second[kXDimension]
  };
}

inline double Dot(const Vector_t &first, const Vector_t &second) {
  return first[kXDimension] * second[kXDimension]
      + first[kYDimension] * second[kYDimension]
      + first[kZDimension] * second[kZDimension];
}

inline double Inner(const Vector_t &vector) {
  return vector[kXDimension] * vector[kXDimension]
      + vector[kYDimension] * vector[kYDimension]
      + vector[kZDimension] * vector[kZDimension];
}

inline Vector_t ElementProduct(const Vector_t &first, const Vector_t &second) {
  return {
      first[kXDimension] * second[kXDimension],
      first[kYDimension] * second[kYDimension],
      first[kZDimension] * second[kZDimension]
  };
}

inline Vector_t ElementDivide(const Vector_t &dividend,
                              const Vector_t &divisor) {
  return {
      dividend[kXDimension] / divisor[kXDimension],
      dividend[kYDimension] / divisor[kYDimension],
      dividend[kZDimension] / divisor[kZDimension]
  };
}

inline double ScalarLength(const Vector_t &vector) {
  return std::sqrt(Inner(vector));
}

inline Vector_t Normalize(const Vector_t &vector) {
  double factor = 1.0 / ScalarLength(vector);
  return vector * factor;
}

inline std::ostream &operator<<(std::ostream &os, const Matrix_t &matrix) {
  os << matrix[kXDimension] << '\n' << matrix[kYDimension] << '\n'
     << matrix[kZDimension];
  return os;
}

inline std::istream &operator>>(std::istream &is, Matrix_t &matrix) {
  is >> matrix[kXDimension] >> matrix[kYDimension] >> matrix[kZDimension];
  return is;
}

inline Matrix_t &operator*=(Matrix_t &lhs, double factor) {
  lhs[kXDimension] *= factor;
  lhs[kYDimension] *= factor;
  lhs[kZDimension] *= factor;
  return lhs;
}

inline Matrix_t &operator/=(Matrix_t &lhs, double divisor) {
  lhs[kXDimension] /= divisor;
  lhs[kYDimension] /= divisor;
  lhs[kZDimension] /= divisor;
  return lhs;
}

inline Matrix_t operator*(const Matrix_t &matrix, double factor) {
  Matrix_t temp(matrix);
  return (temp *= factor);
}

inline Matrix_t operator*(double factor, const Matrix_t &matrix) {
  return operator*(matrix, factor);
}

inline Matrix_t operator/(const Matrix_t &matrix, double divisor) {
  Matrix_t temp(matrix);
  return (temp /= divisor);
}

inline Vector_t operator*(const Vector_t &lhs, const Matrix_t &rhs) {
  return {
      lhs[kXDimension] * rhs[kXDimension][kXDimension]
          + lhs[kYDimension] * rhs[kYDimension][kXDimension]
          + lhs[kZDimension] * rhs[kZDimension][kXDimension],
      lhs[kXDimension] * rhs[kXDimension][kYDimension]
          + lhs[kYDimension] * rhs[kYDimension][kYDimension]
          + lhs[kZDimension] * rhs[kZDimension][kYDimension],
      lhs[kXDimension] * rhs[kXDimension][kZDimension]
          + lhs[kYDimension] * rhs[kYDimension][kZDimension]
          + lhs[kZDimension] * rhs[kZDimension][kZDimension]
  };
}
inline double Determinant(const Matrix_t &input) {
  return (input[kXDimension][kXDimension] * input[kYDimension][kYDimension]
      * input[kZDimension][kZDimension]
      - input[kXDimension][kXDimension] * input[kYDimension][kZDimension]
          * input[kZDimension][kYDimension]
      - input[kXDimension][kYDimension] * input[kYDimension][kXDimension]
          * input[kZDimension][kZDimension]
      + input[kXDimension][kYDimension] * input[kYDimension][kZDimension]
          * input[kZDimension][kXDimension]
      + input[kXDimension][kZDimension] * input[kYDimension][kXDimension]
          * input[kZDimension][kYDimension]
      - input[kXDimension][kZDimension] * input[kYDimension][kYDimension]
          * input[kZDimension][kXDimension]);
}
inline Matrix_t TransposeMatrix(const Matrix_t &input) {
  return {{{input[kXDimension][kXDimension], input[kYDimension][kXDimension],
            input[kZDimension][kXDimension]},
           {input[kXDimension][kYDimension], input[kYDimension][kYDimension],
            input[kZDimension][kYDimension]},
           {input[kXDimension][kZDimension], input[kYDimension][kZDimension],
            input[kZDimension][kZDimension]}}};
}
inline Matrix_t InverseMatrix(const Matrix_t &input) {
  double det = Determinant(input);
  return {
      {
          {
              (input[kYDimension][kYDimension] * input[kZDimension][kZDimension]
                  - input[kYDimension][kZDimension] * input[kZDimension][kYDimension]) / det,
              (input[kXDimension][kZDimension] * input[kZDimension][kYDimension]
                  - input[kXDimension][kYDimension] * input[kZDimension][kZDimension]) / det,
              (input[kXDimension][kYDimension] * input[kYDimension][kZDimension]
                  - input[kXDimension][kZDimension] * input[kYDimension][kYDimension]) / det
          },
          {
              (input[kYDimension][kZDimension] * input[kZDimension][kXDimension]
                  - input[kYDimension][kXDimension] * input[kZDimension][kZDimension]) / det,
              (input[kXDimension][kXDimension] * input[kZDimension][kZDimension]
                  - input[kXDimension][kZDimension] * input[kZDimension][kXDimension]) / det,
              (input[kXDimension][kZDimension] * input[kYDimension][kXDimension]
                  - input[kXDimension][kXDimension] * input[kYDimension][kZDimension]) / det
          },
          {
              (input[kYDimension][kXDimension] * input[kZDimension][kYDimension]
                  - input[kYDimension][kYDimension] * input[kZDimension][kXDimension]) / det,
              (input[kXDimension][kYDimension] * input[kZDimension][kXDimension]
                  - input[kXDimension][kXDimension] * input[kZDimension][kYDimension]) / det,
              (input[kXDimension][kXDimension] * input[kYDimension][kYDimension]
                  - input[kXDimension][kYDimension] * input[kYDimension][kXDimension]) / det
          }
      }
  };
}
#endif //LMC_CFG_INCLUDE_VECTORMATRIX_HPP_
