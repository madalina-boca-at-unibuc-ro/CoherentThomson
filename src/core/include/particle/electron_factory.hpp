#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../io_utils/io_utils.hpp"
#include "../math_utils/math_constants.hpp"
#include "../phys_utils/phys_utils.hpp"
#include "electron.hpp"

namespace Core::Particle {

// Cylinder-beam config parsing lives in IoUtils (shared with
// Simulation::init_simulation_parameters, which needs the same average_px/py/pz).
using CylinderBeamParams = Core::IoUtils::CylinderBeamParams;

// Draws one electron: position uniform over the cylinder cross-section
// (area-preserving in r) and uniform along the cylinder's height, momentum
// Gaussian about params.average_p*. The cylinder's height axis is
// reoriented via beam_axis_rotation to match the beam's own mean momentum
// direction (so the finite-height extent lies along the direction the beam
// actually moves in, rather than always along canonical Oz regardless of
// where average_p* points). Advances gen's state, so a caller generating a
// whole beam must share one generator across calls to get a single
// reproducible sequence.
Electron generate_electron(const CylinderBeamParams& params,
                           const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation, double tau_0, double d_tau,
                           size_t N_tau, std::mt19937& gen);

// Bunched stack, derived from the beam_modulation_* keys. Lengths are in a.u.; z is before beam_center and the beam
// rotation. Sheets start at the uniform z_j = (j - (M-1)/2) d and move to z_j + (zeta / K) sin(K z_j + phase), each
// keeping charge 1/M.
struct ModulationGeometry {
  double q_axial;       // the harmonic-1 translation wavevector along the beam axis, signed
  double q_transverse;  // its size across the beam axis, which equal sheets can't match
  double lambda_1;      // 2 pi / |q_axial|
  double lambda_mod;    // period_ratio * lambda_1
  double K;             // 2 pi / lambda_mod
  double spacing;       // d = lambda_mod / S, before the sheets move
  double length;        // H = M d, before the sheets move
  double spill;         // how far the moved sheets leave +-H/2; 0 if none does
  std::vector<double> sheet_z;
};

// Uses the configured average_p, as the beam doesn't exist yet; it equals the generated mean at zero momentum spread.
ModulationGeometry modulation_geometry(const CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                                       const std::array<double, 3>& detector_normal);

// One line of the log's check: b at harmonic N against the cold-beam value.
struct ModulationHarmonic {
  size_t N;
  double n;                                  // N * period_ratio
  std::complex<double> sheet_sum;            // (1/M) sum_j exp(+i q_N.z_j), as the form factor sees it
  std::optional<std::complex<double>> cold;  // the same from the map's phase integral; integer n, whole periods only
};
std::vector<ModulationHarmonic> modulation_harmonics(const CylinderBeamParams& params,
                                                     const ModulationGeometry& geometry, size_t N_min,
                                                     size_t N_harmonics);

// What the user should know before trusting the stack: a non-whole number of periods, too few sheets for N_max.
std::vector<std::string> modulation_warnings(const CylinderBeamParams& params, size_t N_max);

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config, const double tau_0, const double d_tau,
                                             const size_t N_tau, const Core::Laser::LaserField& laser,
                                             const std::array<double, 3>& detector_normal);

}  // namespace Core::Particle
