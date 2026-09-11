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

// compute the dressed (ponderomotive drift) electron momentum q = p + (mc)^2 <a^2> / (2 k1.p) * k1,
// giving the mass-shell shift m_eff^2 = m^2(1+<a^2>). The (mc)^2 factor is essential: xi is
// dimensionless, so without it the correction term has the wrong units and (since mc = 137.036 in
// atomic units) ends up ~(mc)^2 too small to have any visible effect.
//
// <a^2> is the cycle-averaged normalized-amplitude-squared, NOT xi^2 = a0^2 (the peak amplitude)
// directly -- q is the electron's *drift* momentum (its trajectory averaged over one laser cycle),
// so it must be built from the cycle-averaged field, not its instantaneous peak. For the on-axis
// carrier A_x=A0*a*cos(phi-alpha), A_y=A0*b*cos(phi-beta) (zeta_1=a*e^{i*alpha}, zeta_2=b*e^{i*beta},
// a^2+b^2=1 after create_laser's normalization), <|A|^2> = A0^2*(a^2*<cos^2>+b^2*<cos^2>) = A0^2/2
// exactly, for ANY a,b with a^2+b^2=1 -- i.e. <a^2> = xi^2/2 regardless of polarization state
// (linear, circular, or elliptical), since the cross term between the two orthogonal components
// never appears in |A|^2=A_x^2+A_y^2 and each squared cosine averages to 1/2 independently. This
// was previously computed as xi^2/(2*k1.p) (i.e. using the peak xi^2 in place of <a^2>), which
// (matched against the independent Python cross-check, ~/Dropbox/work/bin/python/
// Superradiant_Thomson's ScreenGeometry.from_parameters) gives m_eff^2=m^2(1+xi^2), not the
// documented m^2(1+xi^2/2) target above -- fixed by using <a^2>=xi^2/2, i.e. a 4, not 2, denominator.
// Independently re-confirmed via a dense_frequency_spectrum=true scan: with <a^2>=xi^2/2 the observed
// spectral peak sits close to omega/omega_1=1 (offset only by the pulse's finite-duration spectral
// width); using xi^2 instead shifts the peak further off omega/omega_1=1, i.e. makes it worse, not
// better -- see CLAUDE.md's "k1, p/q, and n2" bullet.
//
// p is the electron's (average) four-momentum, k1 the incident light-like four-vector (canonical-frame
// omega_laser/c along Oz), xi the laser's peak normalized amplitude (laser.get_a0()).
inline MathUtils::RealFourVector dressed_momentum(const MathUtils::RealFourVector& p,
                                                   const MathUtils::RealFourVector& k1, double xi) {
  double mc = AtomicUnits::m_0 * AtomicUnits::c;
  double a_sq_avg = xi * xi / 2.0;
  return p + mc * mc * a_sq_avg / (2.0 * MathUtils::contract(p, k1)) * k1;
}

}  // namespace Core::PhysUtils
