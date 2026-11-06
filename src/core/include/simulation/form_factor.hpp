#pragma once
#include <complex>
#include <cstddef>
#include <vector>

#include "../detector/detector.hpp"
#include "../laser/laser_field.hpp"
#include "../math_utils/math_utils.hpp"
#include "../particle/electron.hpp"

namespace Core::Simulation {

// Single-wavevector form factor of the initial electron positions (theory/thomson_form_factor_diagnostic.md):
// b(q) = (1/W) sum_j w_j exp(i q.r_j(0)), W = sum_j w_j, with one run-level q (below) and w_j each electron's
// charge weight (1 unless the beam is a lattice stack). A diagnostic only: it does not enter the radiation calculation.
// All wavevectors are four-vectors in atomic units (1/bohr), K = (omega/c, k).
struct FormFactor {
  size_t num_electrons;
  double total_charge_weight;  // W, the charge the simulated electrons stand for, in units of one electron
  // Mean initial momentum: the charge-weighted beam average of the 3-momenta, put back on the mass shell.
  MathUtils::RealFourVector mean_momentum;
  // The same ponderomotively dressed momentum (PhysUtils::dressed_momentum) the radiation frequencies are built from.
  MathUtils::RealFourVector dressed_momentum;
  double omega_s;                 // scattered fundamental, non_linear_Thomson_formula(k_L, dressed, n_s, 1) * c
  MathUtils::RealFourVector k_L;  // incident laser, omega_L/c (1, n_L)
  MathUtils::RealFourVector k_s;  // scattered photon, omega_s/c (1, n_s), n_s = the detector's own normal direction
  // q = k_s - (K_s.P / K_L.P) k_L with the undressed mean P; only its spatial part enters b. Translating an electron
  // by r multiplies its field by exp(-i q.r), the conjugate of b's own phase (theory/dev_notes.md).
  MathUtils::RealFourVector q;
  // k_s - k_L, for comparison: 7e-6 / 3e-5 relative from q head-on at a0 = 0.1 / 0.2, ~2% of the transverse q when
  // tilted.
  MathUtils::RealFourVector q_without_ratio;
  std::complex<double> b;
};

// One pass over the already-generated beam's stored initial positions/momenta (the generator is left untouched).
FormFactor compute_form_factor(const std::vector<Particle::Electron>& electron_beam, const Laser::LaserField& laser,
                               const Detector::Detector_2D& detector);

}  // namespace Core::Simulation
