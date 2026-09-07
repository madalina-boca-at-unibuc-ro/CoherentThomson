#pragma once

#include "../math_utils/math_utils.hpp"
#include "../particle/electron.hpp"
#include "../simulation/simulation.hpp"

namespace Core::Radiation {

// Computes one electron's contribution to the radiated field and adds it into `field` (on top of
// whatever it already holds) at every configured frequency/screen-point pair. `field` holds the
// packed accumulation-time representation (Simulation::PackedFaraday); run_simulation
// reconstructs the true Faraday tensor (Simulation::Faraday) once every electron's and every
// thread's contribution has been summed.
//
// `use_direct_formula` selects which of the two analytically-equivalent closed forms in
// theory/FT_Faraday_tensor-direct_and_simplified_forms.md is evaluated per trajectory point:
// false (default, matches the config key "radiation_formula"="simplified") uses the integration-
// by-parts "simplified" form (an explicit i*omega factor on the long-range term, no acceleration
// needed); true ("direct") uses the direct-FT form instead, which needs the electron's
// 4-acceleration (Particle::Electron::State::acceleration) and has no explicit frequency factor
// outside the phase. Only the *production* path (this function) supports both forms --
// debug/debug_radiation.cpp's diagnostic still only implements the simplified form (see CLAUDE.md).
void compute_radiation(Particle::Electron& electron, const Laser::LaserField& laser,
                       const std::vector<double>& frequencies_list, const Detector::Detector_2D& detector,
                       Simulation::PackedRadiationField& field, bool use_direct_formula);

}  // namespace Core::Radiation
