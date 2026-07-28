#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <type_traits>
#include <utility>

#include "math_constants.hpp"
namespace Core::MathUtils {

template <typename T>
class FourVector;
template <typename T>
class FourTensor;

// The Minkowski Metric Tensor components: g_00 = 1, g_11 = g_22 = g_33 = -1
// Stored as a flat array layout for fast diagonal mappings
inline constexpr std::array<double, 4> g{1.0, -1.0, -1.0, -1.0};

// =========================================================================
// 1. Relativistic FourVector Class (Rank-1 Tensor)
// =========================================================================
template <typename T>
class FourVector {
private:
  // Always stored internally as Contravariant: V^\mu = (V^0, V^1, V^2, V^3)
  std::array<T, 4> data;

public:
  // Constructors
  constexpr FourVector() : data{T(0), T(0), T(0), T(0)} {}
  constexpr FourVector(T v0, T v1, T v2, T v3) : data{v0, v1, v2, v3} {}
  constexpr explicit FourVector(const std::array<T, 4>& raw_data) : data(raw_data) {}

  // Returns natural upper-index components: V^\mu
  constexpr std::array<T, 4> u() const { return data; }

  // Returns lowered-index components on the fly: V_\mu = g_{\mu\nu} V^\nu
  constexpr std::array<T, 4> d() const { return {data[0] * g[0], data[1] * g[1], data[2] * g[2], data[3] * g[3]}; }

  // Subscript operators for direct mutation of contravariant components
  constexpr T& operator[](size_t index) { return data[index]; }
  constexpr const T& operator[](size_t index) const { return data[index]; }

  // copy: create a new FourVector with the same data
  constexpr FourVector<T> copy() const { return FourVector<T>(data); }

  // In-place Vector Addition (V1 += V2)
  constexpr FourVector<T>& operator+=(const FourVector<T>& rhs) noexcept {
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    data[3] += rhs.data[3];
    return *this;
  }

  // In-place Vector Subtraction (V1 -= V2)
  constexpr FourVector<T>& operator-=(const FourVector<T>& rhs) noexcept {
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    data[3] -= rhs.data[3];
    return *this;
  }

  // In-place Scalar Multiplication (V *= scalar)
  // Templated on ScalarT to seamlessly handle multiplication by int, double, or
  // complex
  template <typename ScalarT>
  constexpr FourVector<T>& operator*=(ScalarT scalar) noexcept {
    data[0] *= static_cast<T>(scalar);
    data[1] *= static_cast<T>(scalar);
    data[2] *= static_cast<T>(scalar);
    data[3] *= static_cast<T>(scalar);
    return *this;
  }

  // In-place Scalar Division (V /= scalar)
  // Templated on ScalarT to seamlessly handle division by int, double, or
  // complex
  template <typename ScalarT>
  constexpr FourVector<T>& operator/=(ScalarT scalar) noexcept {
    data[0] /= static_cast<T>(scalar);
    data[1] /= static_cast<T>(scalar);
    data[2] /= static_cast<T>(scalar);
    data[3] /= static_cast<T>(scalar);
    return *this;
  }
};

// Global Out-of-place Vector Scalar Multiplication (V_new = V * scalar)
template <typename T, typename ScalarT>
inline constexpr FourVector<T> operator*(FourVector<T> vector, ScalarT scalar) noexcept {
  vector *= scalar;
  return vector;
}

// Global Out-of-place Vector Scalar Multiplication (V_new = V * scalar)
template <typename T, typename ScalarT>
inline constexpr FourVector<T> operator*(ScalarT scalar, FourVector<T> vector) noexcept {
  vector *= scalar;
  return vector;
}

// Global Out-of-place Vector Scalar Division (V_new = V / scalar)
template <typename T, typename ScalarT>
inline constexpr FourVector<T> operator/(FourVector<T> vector, ScalarT scalar) noexcept {
  vector /= scalar;
  return vector;
}

// Global Out-of-place Binary Operator (V_new = V1 + V2)
template <typename T>
inline constexpr FourVector<T> operator+(FourVector<T> lhs, const FourVector<T>& rhs) noexcept {
  lhs += rhs;  // Efficiently invokes the in-place member operator+= on the local
               // copy
  return lhs;
}

// Global Out-of-place Binary Operator (V_new = V1 - V2)
template <typename T>
inline constexpr FourVector<T> operator-(FourVector<T> lhs, const FourVector<T>& rhs) noexcept {
  lhs -= rhs;  // Efficiently invokes the in-place member operator-= on the local
               // copy
  return lhs;
}

// =========================================================================
// 2. Relativistic FourTensor Class (Rank-2 Tensor Matrix)
// =========================================================================
template <typename T>
class FourTensor {
private:
  // Always stored internally as Fully Contravariant: T^{\mu\nu}
  std::array<std::array<T, 4>, 4> data;

public:
  // Constructors
  constexpr FourTensor() : data{} {}  // Initializes all 16 elements to zero
  constexpr FourTensor(const std::array<std::array<T, 4>, 4>& raw_data) : data(raw_data) {}

  // Both Up: T^{\mu\nu}
  constexpr std::array<std::array<T, 4>, 4> uu() const { return data; }

  // Both Down: T_{\mu\nu} = T^{\alpha\beta} g_{\mu\alpha} g_{\nu\beta}
  constexpr std::array<std::array<T, 4>, 4> dd() const {
    std::array<std::array<T, 4>, 4> result;
    for (size_t mu = 0; mu < 4; ++mu) {
      for (size_t nu = 0; nu < 4; ++nu) {
        result[mu][nu] = data[mu][nu] * g[mu] * g[nu];
      }
    }
    return result;
  }

  // First Up, Second Down: T^\mu_{\ \nu} = T^{\mu\beta} g_{\nu\beta}
  constexpr std::array<std::array<T, 4>, 4> ud() const {
    std::array<std::array<T, 4>, 4> result;
    for (size_t mu = 0; mu < 4; ++mu) {
      for (size_t nu = 0; nu < 4; ++nu) {
        result[mu][nu] = data[mu][nu] * g[nu];
      }
    }
    return result;
  }

  // First Down, Second Up: T_\mu^{\ \nu} = T_{\alpha\nu} g_{\mu\alpha}
  constexpr std::array<std::array<T, 4>, 4> du() const {
    std::array<std::array<T, 4>, 4> result;
    for (size_t mu = 0; mu < 4; ++mu) {
      for (size_t nu = 0; nu < 4; ++nu) {
        result[mu][nu] = data[mu][nu] * g[mu];
      }
    }
    return result;
  }

  // Row/Column mutations for the underlying contravariant matrix data
  constexpr std::array<T, 4>& operator[](size_t index) { return data[index]; }
  constexpr const std::array<T, 4>& operator[](size_t index) const { return data[index]; }

  // In-place Tensor Addition (T1 += T2)
  constexpr FourTensor<T>& operator+=(const FourTensor<T>& rhs) noexcept {
    for (size_t mu = 0; mu < 4; ++mu) {
      for (size_t nu = 0; nu < 4; ++nu) {
        data[mu][nu] += rhs.data[mu][nu];
      }
    }
    return *this;
  }

  // In-place Scalar Multiplication (T *= scalar)
  // Templated on ScalarT to seamlessly handle multiplication by int, double, or
  // complex
  template <typename ScalarT>
  constexpr FourTensor<T>& operator*=(ScalarT scalar) noexcept {
    for (size_t mu = 0; mu < 4; ++mu) {
      for (size_t nu = 0; nu < 4; ++nu) {
        data[mu][nu] *= static_cast<T>(scalar);
      }
    }
    return *this;
  }
};

// Clean Type Aliases for Daily Physics Dev
using Complex = std::complex<double>;
using RealFourVector = FourVector<double>;
using ComplexFourVector = FourVector<std::complex<double>>;
using RealFourTensor = FourTensor<double>;
using ComplexFourTensor = FourTensor<std::complex<double>>;

// A packed antisymmetric bivector: the 6 independent upper-triangle elements of a 4x4
// antisymmetric tensor (diagonal and lower triangle are always zero/mirrored, so never worth
// storing or summing), in the fixed index order (0,1),(0,2),(0,3),(1,2),(1,3),(2,3). Used as a
// cheap accumulation-time stand-in for a ComplexFourTensor known to be antisymmetric -- see
// unpack_bivector below to reconstruct the full tensor once accumulation is finished.
using ComplexBivector = std::array<Complex, 6>;

inline ComplexBivector& operator+=(ComplexBivector& lhs, const ComplexBivector& rhs) noexcept {
  for (size_t i = 0; i < 6; ++i)
    lhs[i] += rhs[i];
  return lhs;
}

// =========================================================================
// 3. Global Math Operations (Dot, Cross, Tensor Contractions)
// =========================================================================

// 3D Spatial Scalar Product (Uses indices 1, 2, 3)
template <typename T>
inline constexpr auto dot3(const FourVector<T>& a, const FourVector<T>& b) {
  auto a_up = a.u();
  auto b_up = b.u();
  if constexpr (std::is_same_v<T, double>) {
    return a_up[1] * b_up[1] + a_up[2] * b_up[2] + a_up[3] * b_up[3];
  } else {
    return std::conj(a_up[1]) * b_up[1] + std::conj(a_up[2]) * b_up[2] + std::conj(a_up[3]) * b_up[3];
  }
}

// 3D Spatial Cross Product
template <typename T>
inline constexpr FourVector<T> cross3(const FourVector<T>& a, const FourVector<T>& b) {
  auto a_up = a.u();
  auto b_up = b.u();
  return FourVector<T>(T(0), a_up[2] * b_up[3] - a_up[3] * b_up[2], a_up[3] * b_up[1] - a_up[1] * b_up[3],
                       a_up[1] * b_up[2] - a_up[2] * b_up[1]);
}

// In place 3D Spatial Cross Product
template <typename T>
inline void cross3_in_place(const FourVector<T>& a, const FourVector<T>& b, FourVector<T>& result) {
  result[0] = T(0);
  result[1] = a[2] * b[3] - a[3] * b[2];
  result[2] = a[3] * b[1] - a[1] * b[3];
  result[3] = a[1] * b[2] - a[2] * b[1];
  return;
}

// =========================================================================
// Overload: FourVector + FourTensor Complex Conjugation
// =========================================================================

// Four-Vector Complex Conjugation
template <typename T>
inline constexpr FourVector<T> conj4(const FourVector<T>& v) {
  if constexpr (std::is_same_v<T, double>) {
    return v;
  } else {
    auto v_up = v.u();
    return FourVector<T>(std::conj(v_up[0]), std::conj(v_up[1]), std::conj(v_up[2]), std::conj(v_up[3]));
  }
}

// Four-Vector Complex Conjugation
template <typename T>
inline constexpr FourTensor<T> conj4(const FourTensor<T>& Tensor) {
  if constexpr (std::is_same_v<T, double>) {
    // A real tensor is its own conjugate! Perfect no-op optimization.
    return Tensor;
  } else {
    std::array<std::array<T, 4>, 4> result;
    auto T_uu = Tensor.uu();

    for (size_t mu = 0; mu < 4; ++mu) {
      for (size_t nu = 0; nu < 4; ++nu) {
        result[mu][nu] = std::conj(T_uu[mu][nu]);
      }
    }
    return FourTensor<T>(result);
  }
}

// =========================================================================
// Overload: FourVector / FourTensor Contraction
// =========================================================================

// Minkowski Four-Product: A^\mu B_\mu
template <typename T>
inline constexpr auto contract(const FourVector<T>& a, const FourVector<T>& b) {
  auto a_up = a.u();
  auto b_low = b.d();
  if constexpr (std::is_same_v<T, double>) {
    return a_up[0] * b_low[0] + a_up[1] * b_low[1] + a_up[2] * b_low[2] + a_up[3] * b_low[3];
  } else {
    return std::conj(a_up[0]) * b_low[0] + std::conj(a_up[1]) * b_low[1] + std::conj(a_up[2]) * b_low[2] +
           std::conj(a_up[3]) * b_low[3];
  }
}

// Tensor Contraction: R^\mu = T^{\mu\nu} v_\nu
template <typename T>
inline constexpr FourVector<T> contract(const FourTensor<T>& Tensor, const FourVector<T>& v) {
  std::array<T, 4> result = {T(0), T(0), T(0), T(0)};
  auto T_uu = Tensor.uu();
  auto v_cov = v.d();  // Grab Covariant v_\nu on the fly to contract correctly

  for (size_t mu = 0; mu < 4; ++mu) {
    result[mu] = T_uu[mu][0] * v_cov[0] + T_uu[mu][1] * v_cov[1] + T_uu[mu][2] * v_cov[2] + T_uu[mu][3] * v_cov[3];
  }
  return FourVector<T>(result);
}

// Tensor Contraction: R^\mu = T^{\mu\nu} v_\nu
template <typename T>
inline constexpr FourVector<T> contract(const FourVector<T>& v, const FourTensor<T>& Tensor) {
  std::array<T, 4> result = {T(0), T(0), T(0), T(0)};
  auto T_uu = Tensor.uu();
  auto v_cov = v.d();  // Grab Covariant v_\nu on the fly to contract correctly

  for (size_t mu = 0; mu < 4; ++mu) {
    result[mu] = T_uu[0][mu] * v_cov[0] + T_uu[1][mu] * v_cov[1] + T_uu[2][mu] * v_cov[2] + T_uu[3][mu] * v_cov[3];
  }
  return FourVector<T>(result);
}

// =========================================================================
// Tensor Transposition: B^{\mu\nu} = A^{\nu\mu}
// Swaps rows and columns of a FourTensor
// =========================================================================
template <typename T>
inline constexpr FourTensor<T> transpose(const FourTensor<T>& Tensor) {
  std::array<std::array<T, 4>, 4> result;
  auto A_uu = Tensor.uu();

  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = 0; nu < 4; ++nu) {
      // Swap the row and column indices
      result[mu][nu] = A_uu[nu][mu];
    }
  }
  return FourTensor<T>(result);
}

// =========================================================================
// Fill in the lower triangle (mu > nu) of an antisymmetric tensor from its upper triangle
// (mu < nu) via T^{\nu\mu} = -T^{\mu\nu}. Assumes the diagonal and lower triangle are still
// at their default-constructed (zero) value on entry -- only the upper-triangle elements
// need to have been written.
// =========================================================================
template <typename T>
inline void mirror_antisymmetric_in_place(FourTensor<T>& tensor) {
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = mu + 1; nu < 4; ++nu) {
      tensor[nu][mu] = -tensor[mu][nu];
    }
  }
}

// Reconstructs the true, full antisymmetric 4x4 tensor T^{mu nu} from its packed upper-triangle
// representation (see ComplexBivector): diagonal is zero, and the lower triangle is filled in as
// -T^{nu mu} via the same index order mirror_antisymmetric_in_place assumes.
inline ComplexFourTensor unpack_bivector(const ComplexBivector& packed) {
  static constexpr std::array<std::pair<size_t, size_t>, 6> index = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  ComplexFourTensor result;
  for (size_t k = 0; k < 6; ++k) {
    auto [alpha, beta] = index[k];
    result[alpha][beta] = packed[k];
    result[beta][alpha] = -packed[k];
  }
  return result;
}

// Type alias for spatial 3D rotations
using RotationMatrix3x3 = std::array<std::array<double, 3>, 3>;

// =========================================================================
// Generate Rotation Matrix: Aligns local +Z (0,0,1) to an arbitrary target
// vector
// =========================================================================
inline constexpr RotationMatrix3x3 rotation_matrix_from_direction(double nx, double ny, double nz) {
  // Normalize the vector to ensure it's a true unit direction
  double norm = std::sqrt(nx * nx + ny * ny + nz * nz);

  constexpr double tolerance = 1e-12;
  if (norm < tolerance) {
    return {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  }

  nx /= norm;
  ny /= norm;
  nz /= norm;

  // Catch edge case: already pointing perfectly along +Z (Identity Matrix)
  if (nz > 0.999999) {
    return {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  }
  // Catch edge case: pointing perfectly along -Z (Flipped Matrix)
  if (nz < -0.999999) {
    return {{{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, -1.0}}};
  }

  // General case using Rodrigues' shortest arc path formula
  double k = 1.0 / (1.0 + nz);
  RotationMatrix3x3 R;
  R[0][0] = 1.0 - nx * nx * k;
  R[0][1] = -nx * ny * k;
  R[0][2] = nx;
  R[1][0] = -nx * ny * k;
  R[1][1] = 1.0 - ny * ny * k;
  R[1][2] = ny;
  R[2][0] = -nx;
  R[2][1] = -ny;
  R[2][2] = nz;

  return R;
}

// =========================================================================
// Apply 3D Rotation to a Spatial vector; rotates an array<double> of size 3
// using a 3x3 rotation matrix
// =========================================================================
inline constexpr std::array<double, 3> rotate3d(const RotationMatrix3x3& R, const std::array<double, 3>& vec) {
  return {R[0][0] * vec[0] + R[0][1] * vec[1] + R[0][2] * vec[2],
          R[1][0] * vec[0] + R[1][1] * vec[1] + R[1][2] * vec[2],
          R[2][0] * vec[0] + R[2][1] * vec[1] + R[2][2] * vec[2]};
}

// =========================================================================
// Compose two 3D rotations: the result applies B first, then A, i.e.
// rotate3d(multiply(A, B), v) == rotate3d(A, rotate3d(B, v))
// =========================================================================
inline constexpr RotationMatrix3x3 multiply(const RotationMatrix3x3& A, const RotationMatrix3x3& B) {
  RotationMatrix3x3 result{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      result[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
    }
  }
  return result;
}

// =========================================================================
// Factorial n!, for small non-negative integer n (e.g. Laguerre-Gauss mode normalization).
// =========================================================================
inline constexpr double factorial(int n) {
  double result = 1.0;
  for (int k = 2; k <= n; ++k) {
    result *= static_cast<double>(k);
  }
  return result;
}

// =========================================================================
// Generalized (associated) Laguerre polynomial L_p^k(x), via the standard three-term recurrence:
// L_0^k = 1, L_1^k = 1 + k - x, (n+1) L_{n+1}^k = (2n+1+k-x) L_n^k - (n+k) L_{n-1}^k.
// Used by the Laguerre-Gauss laser mode's radial profile (p, k = |l| are small non-negative integers
// there, so the recurrence is numerically well-behaved).
// =========================================================================
inline constexpr double generalized_laguerre(int p, int k, double x) {
  if (p == 0) return 1.0;
  double L_prev = 1.0;          // L_0^k
  double L_curr = 1.0 + k - x;  // L_1^k
  for (int n = 1; n < p; ++n) {
    double L_next = ((2 * n + 1 + k - x) * L_curr - (n + k) * L_prev) / (n + 1);
    L_prev = L_curr;
    L_curr = L_next;
  }
  return L_curr;
}

// =========================================================================
// Confluent hypergeometric (Kummer) function 1F1(a; b; z), restricted to a a non-positive integer
// (a = -n, n >= 0). The defining series
//   1F1(a; b; z) = sum_{n=0}^inf [(a)_n / (b)_n] * z^n / n!     ((a)_n, (b)_n: rising factorials)
// then terminates exactly at n = -a (every later term carries a factor (a + n - 1) that is exactly
// zero once n = -a + 1), so this evaluates it as a finite polynomial with no truncation error, via
// the standard term-ratio recurrence t_n = t_{n-1} * (a+n-1)/(b+n-1) * z/n, t_0 = 1.
//
// This is the case needed for the generalized Laguerre polynomials used by the Laguerre-Gauss laser
// mode: L_p^k(x) = binom(p+k, p) * 1F1(-p; k+1; x).
//
// Precondition: a <= 0 and integral. Passing a positive or non-integer a will silently sum only
// the first -a (rounded-down) terms instead of the true (generally infinite, non-terminating) series.
// =========================================================================
inline constexpr double hypergeometric_1F1_neg_int_a(int a, double b, double z) {
  int n_terms = -a;  // number of nonzero terms beyond the leading 1
  double term = 1.0;
  double sum = 1.0;
  for (int n = 1; n <= n_terms; ++n) {
    term *= (static_cast<double>(a + n - 1) / (b + static_cast<double>(n - 1))) * (z / static_cast<double>(n));
    sum += term;
  }
  return sum;
}

// =========================================================================
// d/dz of 1F1(a; b; z), via the standard Kummer derivative identity
//   d/dz 1F1(a; b; z) = (a/b) * 1F1(a+1; b+1; z),
// specialized to the same a <= 0 integer case as hypergeometric_1F1_neg_int_a above (a+1 is then
// still <= 0, so the right-hand side stays within that function's domain). When a = 0, 1F1(0; b; z)
// is the constant 1, so the derivative is 0 -- handled directly rather than by evaluating
// 1F1(1; b+1; z), which would fall outside the negative-integer-only case this header supports.
//
// Precondition: same as hypergeometric_1F1_neg_int_a -- a <= 0 and integral.
// =========================================================================
inline constexpr double hypergeometric_1F1_neg_int_a_derivative(int a, double b, double z) {
  if (a == 0) return 0.0;
  return (static_cast<double>(a) / b) * hypergeometric_1F1_neg_int_a(a + 1, b + 1.0, z);
}

// =========================================================================
// 4x4 Rotation Tensor: Natively stored as Up-Up (T^{\mu\nu})
// =========================================================================
template <typename T = double>
inline constexpr FourTensor<T> rotation_four_tensor_from_direction(double nx, double ny, double nz) {
  // 1. Get the standard 3x3 rotation matrix (which acts as a mixed R^\mu_\nu)
  RotationMatrix3x3 R3 = rotation_matrix_from_direction(nx, ny, nz);

  // 2. Build the 4x4 layout, explicitly multiplying the spatial columns
  //    by g[nu] (which is -1) to raise the second index from R^\mu_\nu to
  //    R^{\mu\nu}.
  std::array<std::array<T, 4>, 4> mat4 = {
      {{T(1.0) * g[0], T(0.0) * g[1], T(0.0) * g[2], T(0.0) * g[3]},
       {T(0.0) * g[0], T(R3[0][0]) * g[1], T(R3[0][1]) * g[2], T(R3[0][2]) * g[3]},
       {T(0.0) * g[0], T(R3[1][0]) * g[1], T(R3[1][1]) * g[2], T(R3[1][2]) * g[3]},
       {T(0.0) * g[0], T(R3[2][0]) * g[1], T(R3[2][1]) * g[2], T(R3[2][2]) * g[3]}}};

  // Because g = {1.0, -1.0, -1.0, -1.0}, this cleanly injects the minus signs
  // into the spatial columns, making the object a true, fully contravariant
  // FourTensor
  return FourTensor<T>(mat4);
};

// =========================================================================
// Inverts a pure-rotation FourTensor built by rotation_four_tensor_from_direction (e.g. a laser's
// rotation_matrix, which maps canonical-frame components -> lab-frame components when used with
// contract()/rotate_tensor() below). A spatial rotation is orthogonal, so its mixed-index form R^mu_nu
// inverts via plain transpose: (R^{-1})^mu_nu = R^nu_mu. R.ud() recovers that ordinary (metric-free)
// mixed matrix (see the comment above rotation_four_tensor_from_direction); the result is re-embedded
// in the same contravariant storage convention as the input so it composes with contract()/
// rotate_tensor() the same way the forward rotation does.
// =========================================================================
template <typename T>
inline constexpr FourTensor<T> inverse_rotation_tensor(const FourTensor<T>& R) {
  auto R_ud = R.ud();
  std::array<std::array<T, 4>, 4> inv_uu{};
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = 0; nu < 4; ++nu) {
      inv_uu[mu][nu] = R_ud[nu][mu] * g[nu];  // transpose, then re-lower nu to restore T^{mu nu} storage
    }
  }
  return FourTensor<T>(inv_uu);
}

// =========================================================================
// Applies the full rank-2 tensor transformation law T'^{alpha beta} = R^alpha_mu R^beta_nu T^{mu nu}
// to a (possibly complex) contravariant tensor, for a rotation FourTensor R in the same convention as
// rotation_four_tensor_from_direction. Pass a laser's own rotation_matrix to rotate canonical -> lab
// (matching contract()'s convention for vectors), or inverse_rotation_tensor(rotation_matrix) to rotate
// lab -> canonical -- e.g. to express a field computed in the lab frame back in the laser's own
// canonical frame (laser along Oz) for plotting.
// =========================================================================
template <typename T, typename U>
inline FourTensor<T> rotate_tensor(const FourTensor<U>& R, const FourTensor<T>& tensor) {
  auto Lambda = R.ud();  // ordinary mixed rotation matrix R^mu_nu, metric factors already cancelled
  auto raw = tensor.uu();

  std::array<std::array<T, 4>, 4> temp{};  // temp^{alpha nu} = R^alpha_mu T^{mu nu}
  for (size_t alpha = 0; alpha < 4; ++alpha) {
    for (size_t nu = 0; nu < 4; ++nu) {
      T sum(0);
      for (size_t mu = 0; mu < 4; ++mu)
        sum += static_cast<T>(Lambda[alpha][mu]) * raw[mu][nu];
      temp[alpha][nu] = sum;
    }
  }

  std::array<std::array<T, 4>, 4> result{};  // T'^{alpha beta} = temp^{alpha nu} R^beta_nu
  for (size_t alpha = 0; alpha < 4; ++alpha) {
    for (size_t beta = 0; beta < 4; ++beta) {
      T sum(0);
      for (size_t nu = 0; nu < 4; ++nu)
        sum += temp[alpha][nu] * static_cast<T>(Lambda[beta][nu]);
      result[alpha][beta] = sum;
    }
  }

  return FourTensor<T>(result);
}

// =========================================================================
// Create a unit light-like 4-vector from an arbitrary input 4-vector by
// normalizing the spatial components and setting the time component to 1.0
// =========================================================================
template <typename T>
inline constexpr FourVector<T> create_unit_light_like_vector(const FourVector<T> vec) {
  FourVector<T> new_vec = vec.copy();  // Create a copy to avoid mutating the original input vector

  double norm = std::sqrt(std::abs(new_vec[1]) * std::abs(new_vec[1]) + std::abs(new_vec[2]) * std::abs(new_vec[2]) +
                          std::abs(new_vec[3]) * std::abs(new_vec[3]));
  if (norm < Core::MathUtils::small_constant) {
    new_vec[0] = 0.0;
    new_vec[1] = 0.0;
    new_vec[2] = 0.0;
    new_vec[3] = 0.0;
  } else {
    new_vec[1] /= norm;
    new_vec[2] /= norm;
    new_vec[3] /= norm;
    new_vec[0] = 1.0;
  }
  return new_vec;
};

// =========================================================================
// Create a unit light-like 4-vector from a direction given by two angles
// =========================================================================
template <typename T>
inline constexpr FourVector<T> create_unit_light_like_vector(const double theta, const double phi) {
  FourVector<T> new_vec{1.0, std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
  return new_vec;
};

template <typename T>
inline double create_unit_light_like_vector_in_place(FourVector<T>& vec) {
  // Returns the original norm of the spatial components before normalization,
  // and mutates the input vector to become a unit light-like vector. If the
  // original norm is below tolerance, the vector is set to zero.
  double norm = std::sqrt(std::abs(vec[1]) * std::abs(vec[1]) + std::abs(vec[2]) * std::abs(vec[2]) +
                          std::abs(vec[3]) * std::abs(vec[3]));
  if (norm < Core::MathUtils::small_constant) {
    vec[0] = 0.0;
    vec[1] = 0.0;
    vec[2] = 0.0;
    vec[3] = 0.0;
  } else {
    double inv_norm = 1.0 / norm;  // Precompute the reciprocal for efficiency
    vec[1] *= inv_norm;
    vec[2] *= inv_norm;
    vec[3] *= inv_norm;
    vec[0] = T(1.0);
  }
  return norm;
};

// Euclidean distance between the spatial (indices 1..3) components of two four-vectors.
template <typename T>
inline constexpr double compute_spatial_separation(const FourVector<T>& a, const FourVector<T>& b) noexcept {
  if constexpr (std::is_same_v<T, double>) {
    double dx = a[1] - b[1];
    double dy = a[2] - b[2];
    double dz = a[3] - b[3];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  } else {
    double dx = std::abs(a[1] - b[1]);
    double dy = std::abs(a[2] - b[2]);
    double dz = std::abs(a[3] - b[3]);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }
};

}  // namespace Core::MathUtils
