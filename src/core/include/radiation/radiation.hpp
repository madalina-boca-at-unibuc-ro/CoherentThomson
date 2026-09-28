#pragma once

#include <string>

#include "../math_utils/math_utils.hpp"
#include "../particle/electron.hpp"
#include "../simulation/simulation.hpp"

namespace Core::Radiation {

// Which closed form of the Fourier-transformed field compute_radiation evaluates, selected by the
// config key "radiation_formula" (parse_radiation_formula):
// - Simplified ("simplified"): exact integration-by-parts form of
//   theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md. Velocity only; fills long_range,
//   short_range, and boundary (the latter only from the two trajectory endpoints).
// - Direct ("direct"): exact direct FT of the Lienard-Wiechert field. Needs the stored
//   4-acceleration; fills long_range and short_range, boundary stays zero.
// - LongDistanceSimplified / LongDistanceDirect ("long_distance_simplified"/"long_distance_direct"):
//   the long-distance approximation of theory/long_distance_direct_simplified_coding_guide.md
//   (fixed reference point r(tau_m), fixed direction n_0, 1/|x_0| amplitude, linearized phase, no
//   1/R^2 term), evaluated by the two-transverse-projection method. short_range stays zero; the
//   simplified variant puts its bulk term in long_range and its endpoint term in boundary, the direct
//   variant everything in long_range. Their long_range/boundary split uses u_perp, so it differs from
//   the exact Simplified split; only the sum is comparable.
// - LongDistanceSimplifiedApprox / LongDistanceDirectApprox ("long_distance_simplified_approx"/
//   "long_distance_direct_approx"): the same formulas and outputs as the two above, but the phase factor
//   e^{ik Phi} is advanced along tau by a small-angle polynomial instead of a sin/cos per (screen point, tau),
//   re-anchored exactly every few dozen steps. Several times faster; agrees with the unstepped variants
//   to ~1e-10 relative.
enum class RadiationFormula {
  Simplified,
  Direct,
  LongDistanceSimplified,
  LongDistanceDirect,
  LongDistanceSimplifiedApprox,
  LongDistanceDirectApprox
};

// Throws std::runtime_error on an unrecognized name, so a typo can't silently select another formula.
RadiationFormula parse_radiation_formula(const std::string& name);

// Computes one electron's contribution to the radiated field and adds it into `field` (on top of
// whatever it already holds) at every configured frequency/screen-point pair. `field` holds the
// packed accumulation-time representation (Simulation::PackedFaraday); run_simulation
// reconstructs the true Faraday tensor (Simulation::Faraday) once every electron's and every
// thread's contribution has been summed.
void compute_radiation(Particle::Electron& electron, const Laser::LaserField& laser,
                       const std::vector<double>& frequencies_list, const Detector::Detector_2D& detector,
                       Simulation::PackedRadiationField& field, RadiationFormula formula);

}  // namespace Core::Radiation
