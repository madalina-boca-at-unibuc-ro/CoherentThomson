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
// (area-preserving in r) and uniform along the cylinder's height, plus the Gaussian wings, then moved by the
// modulation if it is on; momentum Gaussian about params.average_p*. The cylinder's height axis is
// reoriented via beam_axis_rotation to match the beam's own mean momentum
// direction (so the finite-height extent lies along the direction the beam
// actually moves in, rather than always along canonical Oz regardless of
// where average_p* points). Advances gen's state, so a caller generating a
// whole beam must share one generator across calls to get a single
// reproducible sequence. params must be resolved (resolve_beam_lengths).
Electron generate_electron(const CylinderBeamParams& params,
                           const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation, double tau_0, double d_tau,
                           size_t N_tau, std::mt19937& gen);

// The harmonic-1 translation wavevector q along the detector normal (the q of Simulation::compute_form_factor, from
// the configured average_p, undressed), split along and across the beam axis. Uses the configured average_p, as the
// beam doesn't exist yet; it equals the generated mean at zero momentum spread.
struct AxialTranslation {
  double q_axial;       // signed, along the beam axis
  double q_transverse;  // its size across the beam axis, which a modulation along the axis can't match
  bool degenerate;      // q (numerically) perpendicular to the axis: lambda_1 is undefined
  double lambda_1;      // 2 pi / |q_axial|; +infinity when degenerate
};
AxialTranslation axial_translation(const CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                                   const std::array<double, 3>& detector_normal);

// Converts the lengths given in lambda_1 to a.u. and sets params.lambda_1 (when q isn't perpendicular to the axis).
// Throws if a length is in lambda_1 and lambda_1 is undefined. Does nothing to lengths already resolved.
void resolve_beam_lengths(CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                          const std::array<double, 3>& detector_normal);

// The modulation's map, z -> z + (zeta / K) sin(K z + phase), K = 2 pi / period. z is along the beam axis from the
// beam centre, before beam_center and the beam rotation.
double modulated_z(const CylinderBeamParams::Modulation& modulation, double z);

// The width a lattice stack's sheets fill: H, plus sigma_z_cut sigma_z on each side when sigma_z > 0.
double stack_span(const CylinderBeamParams& params);

// A lattice stack's sheet positions: cell-centred in the span L, z_m = (m - (M-1)/2) L/M, moved by the modulation if it
// is on.
std::vector<double> sheet_positions(const CylinderBeamParams& params);

// A lattice stack's sheet charges, summing to 1: exactly 1/M each when sigma_z = 0; otherwise the smoothed density
// (top-hat of width H convolved with N(0, sigma_z^2)) at each unmodulated z_m, normalized. The modulation moves a sheet
// with its charge.
std::vector<double> sheet_weights(const CylinderBeamParams& params);

// Each electron's initial position along the beam axis, from the beam centre: the z the generator drew (and the
// modulation moved), recovered from the rotated position.
std::vector<double> axial_positions(const std::vector<Electron>& beam, const CylinderBeamParams& params);

// The modulation's derived numbers, from resolved params. Lengths are in a.u.
struct ModulationGeometry {
  AxialTranslation q;
  double lambda_mod;  // beam_modulation_period
  double ratio;       // lambda_mod / lambda_1; 0 when lambda_1 is undefined
  double K;           // 2 pi / lambda_mod
  double periods;     // H / lambda_mod
  double span;        // lattice: stack_span, L; random: H
  double spill;       // lattice: how far the moved sheets leave +-L/2 (0 if none does); random: 0
};
ModulationGeometry modulation_geometry(const CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                                       const std::array<double, 3>& detector_normal);

// One line of the log's check: b at harmonic N against the cold-beam value.
struct ModulationHarmonic {
  size_t N;
  double n;                                  // N * ratio
  std::complex<double> sum;                  // (1/count) sum_j exp(+i q_N z_j), as the form factor sees it
  std::optional<std::complex<double>> cold;  // the map's phase integral: integer n, whole periods, no wings only
};
// axial_z: the sheets of a lattice stack, or every electron of a random beam. weights: their charges, summing to 1
// (sheet_weights); empty for equal charges, the sum then being divided by the count.
std::vector<ModulationHarmonic> modulation_harmonics(const CylinderBeamParams& params,
                                                     const ModulationGeometry& geometry,
                                                     const std::vector<double>& axial_z,
                                                     const std::vector<double>& weights, size_t N_min,
                                                     size_t N_harmonics);

// One line of the log's check for an unmodulated Gaussian stack: the weighted sheet sum at harmonic N against the
// continuous density's sinc(q_N H / 2) exp(-q_N^2 sigma_z^2 / 2), q_N = N |q_axial|.
struct StackHarmonic {
  size_t N;
  std::complex<double> sum;  // sum_m w_m exp(+i q_N z_m), signed q_axial as the form factor sees it
  double smooth;             // the continuous density's value (real: the density is even)
};
std::vector<StackHarmonic> stack_harmonics(const CylinderBeamParams& params, const AxialTranslation& q, size_t N_min,
                                           size_t N_harmonics);

// What the user should know about a Gaussian stack's sampling (lattice, M > 1, sigma_z > 0): sheets too far apart for
// the Gaussian edges at N_max, and a tail cut that leaves much density at the span's ends.
std::vector<std::string> stack_warnings(const CylinderBeamParams& params, const AxialTranslation& q, size_t N_max);

// What the user should know before trusting the modulation: a non-whole number of periods, too few sheets per period
// (counted over the span L) for N_max.
std::vector<std::string> modulation_warnings(const CylinderBeamParams& params, const ModulationGeometry& geometry,
                                             size_t N_max);

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config, const double tau_0, const double d_tau,
                                             const size_t N_tau, const Core::Laser::LaserField& laser,
                                             const std::array<double, 3>& detector_normal);

}  // namespace Core::Particle
