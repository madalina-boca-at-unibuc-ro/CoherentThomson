#pragma once
#include "../math_utils/math_constants.hpp"
#include "../math_utils/math_utils.hpp"

namespace Core::PhysUtils::AtomicUnits {
// Fundamental constants in atomic units
inline constexpr double m_0 = 1.0;
inline constexpr double e_0 = 1.0;                                    // the electron charge in modulus
inline constexpr double q_0 = -1.0;                                   // the electron charge
inline constexpr double c = 137.036;                                  // the speed of light
inline constexpr double hbar = 1.0;                                   // the reduced planck constant
inline constexpr double epsilon_0 = 1.0 / (4 * Core::MathUtils::pi);  // the vacuum permittivity
inline constexpr double mu_0 = 1.0 / (epsilon_0 * c * c);             // the vacuum permeability

}  // namespace Core::PhysUtils::AtomicUnits

namespace Core::PhysUtils {

inline double non_linear_Thomson_formula(const MathUtils::RealFourVector& k1, const MathUtils::RealFourVector& p1,
                                         const MathUtils::RealFourVector& n2, size_t N) {
  double omega1 = k1[0] * AtomicUnits::c;
  MathUtils::RealFourVector n1 = k1 / k1[0];
  double omega2 = omega1 * N * MathUtils::contract(p1, n1) / MathUtils::contract(p1, n2);
  return omega2 / AtomicUnits::c;  // Here the frequencies are returned divided by the speed of light, since they appear
                                   // as such in the formula of the FT. TO BE DOCUMENTED IN CLAUDE.md ; the agreement
                                   // with the implementation should be checked
};

}  // namespace Core::PhysUtils
