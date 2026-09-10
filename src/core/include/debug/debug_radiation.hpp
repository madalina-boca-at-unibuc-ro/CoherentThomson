#pragma once

#include <string>

#include "../math_utils/math_utils.hpp"
#include "../particle/electron.hpp"

namespace Core::Debug {

// Debug-only diagnostic, deliberately kept independent of Radiation::compute_radiation: for one
// electron's already-computed trajectory (electron.compute_trajectory() must have been called
// beforehand -- this function only reads electron.get_trajectory(), it never integrates one) and
// one screen point, at a single wavenumber `k` (= omega/c, same convention as
// Radiation::compute_radiation's frequencies_list -- see radiation.cpp's "Here frequencies are
// actually divided by c" comment -- not raw omega), writes one row per trajectory point (indexed by
// tau) of the per-tau long-range/short-range/boundary contributions to the 6 independent Faraday
// bivector components -- i.e. the raw terms Radiation::compute_radiation sums over tau, rather than
// only their final sum -- to `filepath`. The boundary column is a per-*row* value (nonzero only at
// i_tau==0/i_tau==N_tau-1), not itself a sum over tau -- see radiation.cpp's boundary_prefactor doc
// comment and theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Form 2's boundary term F_b"
// section; this function always evaluates the simplified formula's math only (see CLAUDE.md's
// "Debug mode" section), so there is no direct-formula equivalent to keep in sync here.
//
// The underlying physics (n_R0/u geometry, retarded phase, amp_long/amp_short/amp_boundary) is
// intentionally reimplemented here rather than shared with radiation.cpp, so this diagnostic path
// can never be affected by, or accidentally affect, the production accumulation path.
void export_radiation_integrand(const Particle::Electron& electron, const MathUtils::RealFourVector& detector_point,
                                double k, const std::string& filepath);

// Companion to export_radiation_integrand above: writes the one phase factor common to every
// long-range/short-range/boundary bivector component at every (tau, screen point) -- the geometric
// bivector_element(n_R0, u, alpha, beta) terms vary per component and are handled there, but
// exp(i * phase_base * k), phase_base = x[0] + R, multiplies all of them identically -- so it's
// exported separately rather than repeated 18 times (6 components x long/short/boundary) per row of
// export_radiation_integrand's file. Same preconditions and `k` convention as
// export_radiation_integrand above (electron.get_trajectory() must already be populated; k =
// omega/c, not raw omega).
void export_radiation_phase(const Particle::Electron& electron, const MathUtils::RealFourVector& detector_point,
                            double k, const std::string& filepath);

}  // namespace Core::Debug
