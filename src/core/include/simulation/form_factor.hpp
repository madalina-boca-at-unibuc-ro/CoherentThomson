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
// b(q) = (1/N) sum_j exp(i q.r_j(0)), with one run-level q = k_s - k_L. A diagnostic only: it does not enter the
// radiation calculation. All wavevectors are four-vectors in atomic units (1/bohr), K = (omega/c, k).
struct FormFactor {
  size_t num_electrons;
  // Mean initial momentum: the beam average of the 3-momenta, put back on the mass shell.
  MathUtils::RealFourVector mean_momentum;
  // The same ponderomotively dressed momentum (PhysUtils::dressed_momentum) the radiation frequencies are built from.
  MathUtils::RealFourVector dressed_momentum;
  double omega_s;                 // scattered fundamental, non_linear_Thomson_formula(k_L, dressed, n_s, 1) * c
  MathUtils::RealFourVector k_L;  // incident laser, omega_L/c (1, n_L)
  MathUtils::RealFourVector k_s;  // scattered photon, omega_s/c (1, n_s), n_s = the detector's own normal direction
  // q = k_s - k_L; only its spatial part enters b. Whether the k_L term should instead carry the factor
  // (K_s.P)/(K_L.P) once omega_s is dressed is an open question (see CLAUDE.md).
  MathUtils::RealFourVector q;
  std::complex<double> b;
};

// One pass over the already-generated beam's stored initial positions/momenta (the generator is left untouched).
FormFactor compute_form_factor(const std::vector<Particle::Electron>& electron_beam, const Laser::LaserField& laser,
                               const Detector::Detector_2D& detector);

}  // namespace Core::Simulation
