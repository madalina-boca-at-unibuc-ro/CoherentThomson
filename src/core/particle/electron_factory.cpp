#include "../include/particle/electron_factory.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

namespace Core::Particle {

namespace {

// Places one electron at (x, y, z) in the cylinder's canonical frame and draws its momentum: translation by
// beam_center, rotation onto the mean momentum, then Gaussian momenta about average_p (advances gen).
Electron place_electron(const CylinderBeamParams& params, const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation,
                        double x, double y, double z, double tau_0, double d_tau, size_t N_tau, std::mt19937& gen,
                        double charge_weight = 1.0) {
  // Apply the beam's global translation before the rotation below, so the offset itself is expressed in the
  // same canonical frame as beam_cylinder_radius/height and gets carried along with the beam through it (rather
  // than being a fixed shift applied afterwards).
  x += params.center_x;
  y += params.center_y;
  z += params.center_z;

  // Reorient the cylinder so its height axis matches the beam's own mean momentum direction.
  std::array<double, 3> pos_local = MathUtils::rotate3d(beam_axis_rotation, {x, y, z});

  // Create the position 4-vector: (t=tau_0, x, y, z)
  MathUtils::RealFourVector position(PhysUtils::AtomicUnits::c * tau_0, pos_local[0], pos_local[1], pos_local[2]);

  double px = params.average_px + params.sigma_px * std::normal_distribution<double>(0.0, 1.0)(gen);
  double py = params.average_py + params.sigma_py * std::normal_distribution<double>(0.0, 1.0)(gen);
  double pz = params.average_pz + params.sigma_pz * std::normal_distribution<double>(0.0, 1.0)(gen);
  double mc = PhysUtils::AtomicUnits::m_0 * PhysUtils::AtomicUnits::c;
  double p0 = std::sqrt(px * px + py * py + pz * pz + mc * mc);  // Energy component from the relativistic
                                                                 // energy-momentum relation
  MathUtils::RealFourVector momentum(p0, px, py, pz);            // Initial momentum 4-vector

  return Electron(position, momentum, tau_0, d_tau, N_tau, charge_weight);
}

}  // namespace

Electron generate_electron(const CylinderBeamParams& params,
                           const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation, double tau_0, double d_tau,
                           size_t N_tau, std::mt19937& gen) {
  if (!params.lengths_resolved) {
    throw std::logic_error("generate_electron: beam lengths in lambda_1 not resolved (resolve_beam_lengths)");
  }
  std::uniform_real_distribution<double> dis_uniform(0.0, 1.0);
  std::uniform_real_distribution<double> dis_z(-params.height / 2.0, params.height / 2.0);

  // Area-preserving radial distribution transformation
  double r = params.radius * std::sqrt(dis_uniform(gen));
  double phi = dis_uniform(gen) * 2.0 * MathUtils::pi;
  double z = dis_z(gen);
  // Gaussian-smoothed cylinder (theory/gaussian_smoothed_electron_cylinder.md): a uniform z on [-H/2, H/2] plus an
  // independent N(0, sigma_z^2) displacement, giving longitudinal form factor sinc(qH/2) * exp(-q^2 sigma_z^2 / 2).
  // Not clipped, so the tails extend past the cylinder ends. The draw is skipped at sigma_z = 0 so the RNG stream,
  // and hence the whole beam, is identical to the sharp-edged cylinder's.
  if (params.sigma_z > 0.0) {
    z += std::normal_distribution<double>(0.0, params.sigma_z)(gen);
  }
  // The modulation moves the electron by its own position, wings included. It draws nothing, so the RNG stream is
  // the unmodulated beam's.
  if (params.modulation.active) {
    z = modulated_z(params.modulation, z);
  }

  // Convert cylindrical coordinates to Cartesian coordinates (cylinder height axis along canonical Oz).
  return place_electron(params, beam_axis_rotation, r * std::cos(phi), r * std::sin(phi), z, tau_0, d_tau, N_tau, gen);
}

AxialTranslation axial_translation(const CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                                   const std::array<double, 3>& detector_normal) {
  const double c = PhysUtils::AtomicUnits::c;
  const double mc = PhysUtils::AtomicUnits::m_0 * c;
  const double px = params.average_px, py = params.average_py, pz = params.average_pz;
  MathUtils::RealFourVector P(std::sqrt(px * px + py * py + pz * pz + mc * mc), px, py, pz);

  // Same k_s and q as Simulation::compute_form_factor, at harmonic 1 along the detector normal.
  MathUtils::RealFourVector k_L = laser.get_unity_n() * laser.get_omega() / c;
  MathUtils::RealFourVector dressed = PhysUtils::dressed_momentum(P, k_L, laser.get_a0());
  MathUtils::RealFourVector n_s(1.0, detector_normal[0], detector_normal[1], detector_normal[2]);
  MathUtils::RealFourVector k_s = n_s * PhysUtils::non_linear_Thomson_formula(k_L, dressed, n_s, 1);
  MathUtils::RealFourVector q = PhysUtils::translation_wavevector(P, k_L, k_s);

  const std::array<double, 3> axis =
      MathUtils::rotate3d(MathUtils::rotation_matrix_from_direction(px, py, pz), {0.0, 0.0, 1.0});
  AxialTranslation result;
  result.q_axial = q[1] * axis[0] + q[2] * axis[1] + q[3] * axis[2];
  const double q_squared = q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  result.q_transverse = std::sqrt(std::max(0.0, q_squared - result.q_axial * result.q_axial));
  result.degenerate = std::abs(result.q_axial) < 1e-6 * std::sqrt(q_squared);
  result.lambda_1 = result.degenerate ? std::numeric_limits<double>::infinity()
                                      : 2.0 * MathUtils::pi / std::abs(result.q_axial);
  return result;
}

void resolve_beam_lengths(CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                          const std::array<double, 3>& detector_normal) {
  const AxialTranslation q = axial_translation(params, laser, detector_normal);
  params.lambda_1 = q.degenerate ? 0.0 : q.lambda_1;
  if (params.lengths_resolved) {
    return;
  }
  if (q.degenerate) {
    throw std::invalid_argument("unit lambda_1: q is perpendicular to the beam axis, so lambda_1 = 2 pi / |q_axial| "
                                "is undefined");
  }
  if (params.height_in_lambda_1) {
    params.height *= q.lambda_1;
  }
  if (params.sigma_z_in_lambda_1) {
    params.sigma_z *= q.lambda_1;
  }
  if (params.modulation.active && params.modulation.period_in_lambda_1) {
    params.modulation.period *= q.lambda_1;
  }
  params.lengths_resolved = true;
}

double modulated_z(const CylinderBeamParams::Modulation& modulation, double z) {
  const double K = 2.0 * MathUtils::pi / modulation.period;
  return z + modulation.zeta / K * std::sin(K * z + modulation.phase);
}

double stack_span(const CylinderBeamParams& params) {
  return params.sigma_z > 0.0 ? params.height + 2.0 * params.sigma_z_cut * params.sigma_z : params.height;
}

namespace {

std::vector<double> unmodulated_sheet_positions(const CylinderBeamParams& params) {
  const size_t sheets = params.stack_sheets;
  const double spacing = stack_span(params) / static_cast<double>(sheets);
  std::vector<double> sheet_z;
  sheet_z.reserve(sheets);
  for (size_t m = 0; m < sheets; ++m) {
    sheet_z.push_back((static_cast<double>(m) - 0.5 * static_cast<double>(sheets - 1)) * spacing);
  }
  return sheet_z;
}

// The top-hat of width H convolved with N(0, sigma^2), up to a constant factor: a Gaussian at H = 0. Written with
// erfc of |z| so the tails don't come from the difference of two erf values near 1.
double smoothed_density(double z, double height, double sigma) {
  const double a = std::abs(z);
  if (height == 0.0) {
    return std::exp(-a * a / (2.0 * sigma * sigma));
  }
  const double s = std::sqrt(2.0) * sigma;
  return 0.5 * (std::erfc((a - height / 2.0) / s) - std::erfc((a + height / 2.0) / s));
}

}  // namespace

std::vector<double> sheet_positions(const CylinderBeamParams& params) {
  std::vector<double> sheet_z = unmodulated_sheet_positions(params);
  if (params.modulation.active) {
    for (double& z : sheet_z) {
      z = modulated_z(params.modulation, z);
    }
  }
  return sheet_z;
}

std::vector<double> sheet_weights(const CylinderBeamParams& params) {
  const size_t sheets = params.stack_sheets;
  if (params.sigma_z == 0.0) {
    return std::vector<double>(sheets, 1.0 / static_cast<double>(sheets));
  }
  std::vector<double> weights;
  weights.reserve(sheets);
  double total = 0.0;
  for (double z : unmodulated_sheet_positions(params)) {
    weights.push_back(smoothed_density(z, params.height, params.sigma_z));
    total += weights.back();
  }
  for (double& w : weights) {
    w /= total;
  }
  return weights;
}

std::vector<double> axial_positions(const std::vector<Electron>& beam, const CylinderBeamParams& params) {
  const std::array<double, 3> axis = MathUtils::rotate3d(
      MathUtils::rotation_matrix_from_direction(params.average_px, params.average_py, params.average_pz),
      {0.0, 0.0, 1.0});
  std::vector<double> z;
  z.reserve(beam.size());
  for (const Electron& electron : beam) {
    const MathUtils::RealFourVector& x = electron.get_initial_position();
    z.push_back(x[1] * axis[0] + x[2] * axis[1] + x[3] * axis[2] - params.center_z);
  }
  return z;
}

ModulationGeometry modulation_geometry(const CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                                       const std::array<double, 3>& detector_normal) {
  if (!params.lengths_resolved) {
    throw std::logic_error("modulation_geometry: beam lengths in lambda_1 not resolved (resolve_beam_lengths)");
  }
  ModulationGeometry geometry;
  geometry.q = axial_translation(params, laser, detector_normal);
  geometry.lambda_mod = params.modulation.period;
  geometry.ratio = geometry.q.degenerate ? 0.0 : geometry.lambda_mod / geometry.q.lambda_1;
  geometry.K = 2.0 * MathUtils::pi / geometry.lambda_mod;
  geometry.periods = params.height / geometry.lambda_mod;
  geometry.span = params.polar_lattice ? stack_span(params) : params.height;
  geometry.spill = 0.0;
  if (params.polar_lattice) {
    for (double z : sheet_positions(params)) {
      geometry.spill = std::max(geometry.spill, std::abs(z) - geometry.span / 2.0);
    }
  }
  return geometry;
}

namespace {

// x if it is within a relative 1e-9 of a whole number (at least 1 when positive_only), nullopt otherwise.
std::optional<double> nearly_whole(double x, bool positive_only) {
  const double whole = std::round(x);
  if (std::abs(x - whole) > 1e-9 * std::max(1.0, std::abs(x)) || (positive_only && whole < 1.0)) {
    return std::nullopt;
  }
  return whole;
}

}  // namespace

std::vector<ModulationHarmonic> modulation_harmonics(const CylinderBeamParams& params,
                                                     const ModulationGeometry& geometry,
                                                     const std::vector<double>& axial_z,
                                                     const std::vector<double>& weights, size_t N_min,
                                                     size_t N_harmonics) {
  if (!weights.empty() && weights.size() != axial_z.size()) {
    throw std::logic_error("modulation_harmonics: one weight per position, or none");
  }
  const CylinderBeamParams::Modulation& modulation = params.modulation;
  const double two_pi = 2.0 * MathUtils::pi;
  const double sign = geometry.q.q_axial > 0.0 ? 1.0 : -1.0;
  // The cold-beam value is the map's average over whole periods of a uniform density: no Gaussian wings.
  const bool whole_periods = nearly_whole(geometry.periods, true).has_value() && params.sigma_z == 0.0;
  std::vector<ModulationHarmonic> harmonics;
  for (size_t N = N_min; N < N_min + N_harmonics; ++N) {
    const double n = static_cast<double>(N) * geometry.ratio;
    std::complex<double> sum(0.0, 0.0);
    for (size_t j = 0; j < axial_z.size(); ++j) {
      sum += std::polar(weights.empty() ? 1.0 : weights[j], sign * n * geometry.K * axial_z[j]);
    }
    if (weights.empty()) {
      sum /= static_cast<double>(axial_z.size());
    }
    ModulationHarmonic harmonic{N, n, sum, std::nullopt};
    // The midpoint rule on this periodic integrand converges fast once it has more points than n (1 + zeta), so it
    // stands in for (-1)^n J_n(n zeta) without needing a Bessel function from the standard library.
    const std::optional<double> whole_n = nearly_whole(n, true);
    if (whole_periods && whole_n) {
      const double n_int = *whole_n;
      const size_t points = 8 * static_cast<size_t>(std::ceil(n_int * (1.0 + modulation.zeta)) + 64.0);
      std::complex<double> cold(0.0, 0.0);
      for (size_t i = 0; i < points; ++i) {
        const double phi_0 = (static_cast<double>(i) + 0.5) * two_pi / static_cast<double>(points);
        const double phi = phi_0 + modulation.zeta * std::sin(phi_0);
        cold += std::polar(1.0, sign * n_int * (phi - modulation.phase));
      }
      harmonic.cold = cold / static_cast<double>(points);
    }
    harmonics.push_back(harmonic);
  }
  return harmonics;
}

std::vector<std::string> modulation_warnings(const CylinderBeamParams& params, const ModulationGeometry& geometry,
                                             size_t N_max) {
  std::vector<std::string> warnings;
  if (geometry.q.degenerate) {
    warnings.push_back("q is perpendicular to the beam axis: a modulation along the axis can't match it, and there is "
                       "no lambda_1 to compare the period with");
    return warnings;
  }
  if (params.height > 0.0 && !nearly_whole(geometry.periods, true)) {
    std::ostringstream text;
    text << "H / Lambda_mod = " << geometry.periods
         << ", not a whole number of periods: b doesn't reduce to the cold-beam value";
    warnings.push_back(text.str());
  }
  if (params.polar_lattice && params.stack_sheets > 1) {
    // Aliasing of the sheet sum (theory/lattice-stack-implemented.md, "Sheets per period"): S equally spaced moved
    // sheets add the aliases J_{mS-n}(n zeta) to (-1)^n J_n(n zeta), negligible once S - n clears n zeta by the Bessel
    // transition width ~n^(1/3). S >= n (1 + zeta) + 5 n^(1/3) keeps the relative error of b_n below ~1e-4
    // (sheets_per_period_aliasing.py). The sheets are equally spaced over the span L, Gaussian tails included.
    const double sheets_per_period = static_cast<double>(params.stack_sheets) * geometry.lambda_mod / geometry.span;
    const double n_max = static_cast<double>(N_max) * geometry.ratio;
    const double recommended = n_max * (1.0 + params.modulation.zeta) + 5.0 * std::cbrt(n_max);
    if (sheets_per_period < recommended) {
      std::ostringstream text;
      text << sheets_per_period << " sheets per period is below n (1 + zeta) + 5 n^(1/3) = " << recommended
           << " at n = N_max Lambda_mod / Lambda_1 = " << n_max
           << ": aliasing can make b wrong at the highest harmonics; see the log's sum against the cold-beam value";
      warnings.push_back(text.str());
    }
  }
  return warnings;
}

std::vector<StackHarmonic> stack_harmonics(const CylinderBeamParams& params, const AxialTranslation& q, size_t N_min,
                                           size_t N_harmonics) {
  const std::vector<double> sheet_z = sheet_positions(params);
  const std::vector<double> weights = sheet_weights(params);
  std::vector<StackHarmonic> harmonics;
  for (size_t N = N_min; N < N_min + N_harmonics; ++N) {
    const double q_N = static_cast<double>(N) * q.q_axial;
    std::complex<double> sum(0.0, 0.0);
    for (size_t m = 0; m < sheet_z.size(); ++m) {
      sum += std::polar(weights[m], q_N * sheet_z[m]);
    }
    const double x = q_N * params.height / 2.0;
    const double sinc = x == 0.0 ? 1.0 : std::sin(x) / x;
    harmonics.push_back({N, sum, sinc * std::exp(-q_N * q_N * params.sigma_z * params.sigma_z / 2.0)});
  }
  return harmonics;
}

std::vector<std::string> stack_warnings(const CylinderBeamParams& params, const AxialTranslation& q, size_t N_max) {
  std::vector<std::string> warnings;
  if (!params.polar_lattice || params.stack_sheets < 2 || params.sigma_z == 0.0 || q.degenerate) {
    return warnings;
  }
  // Equally spaced samples of the smoothed density add replicas of its transform shifted by 2 pi / d. The nearest one
  // to q_N falls off as exp(-((2 pi / d - q_N) sigma_z)^2 / 2), so (2 pi / d - q_N) sigma_z >= 6 keeps it near e^-18
  // (theory/lattice-stack-implemented.md, "Gaussian stack"). The sinc of a wide H only lowers it further.
  const double span = stack_span(params);
  const double spacing = span / static_cast<double>(params.stack_sheets);
  const double q_max = static_cast<double>(N_max) * std::abs(q.q_axial);
  const double margin = (2.0 * MathUtils::pi / spacing - q_max) * params.sigma_z;
  if (margin < 6.0) {
    const double needed = std::ceil(span * (q_max + 6.0 / params.sigma_z) / (2.0 * MathUtils::pi));
    std::ostringstream text;
    text << "sheet spacing L/M = " << spacing / params.sigma_z
         << " sigma_z is too coarse for the Gaussian edges at N = " << N_max << ": (2 pi / d - q_N) sigma_z = " << margin << " < 6; beam_stack_sheets >= " << needed
         << " keeps the replica below ~1e-8";
    warnings.push_back(text.str());
  }
  if (params.sigma_z_cut < 5.0) {
    std::ostringstream text;
    text << "beam_stack_sigma_z_cut " << params.sigma_z_cut << " cuts the Gaussian where its density is still "
         << std::exp(-params.sigma_z_cut * params.sigma_z_cut / 2.0)
         << " of the peak: the span's ends add their own ripple to b at about that level";
    warnings.push_back(text.str());
  }
  return warnings;
}

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config, const double tau_0, const double d_tau,
                                             const size_t N_tau, const Core::Laser::LaserField& laser,
                                             const std::array<double, 3>& detector_normal) {
  CylinderBeamParams params = IoUtils::parse_cylinder_beam_params(config);
  resolve_beam_lengths(params, laser, detector_normal);
  for (const std::string& note : params.notes) {
    std::cerr << "Warning: " << note << "\n";
  }
  const size_t N_max = IoUtils::get_number_of_harmonics_min(config) + IoUtils::get_number_of_harmonics(config) - 1;
  if (params.modulation.active) {
    for (const std::string& warning :
         modulation_warnings(params, modulation_geometry(params, laser, detector_normal), N_max)) {
      std::cerr << "Warning: " << warning << "\n";
    }
  }
  for (const std::string& warning :
       stack_warnings(params, axial_translation(params, laser, detector_normal), N_max)) {
    std::cerr << "Warning: " << warning << "\n";
  }

  // Aligns the cylinder's height axis with the beam's own mean momentum direction (average_px/py/pz), so the
  // finite-height extent lies along the direction the beam actually moves in rather than always along canonical
  // Oz. Falls back to the identity if the mean momentum is (numerically) zero.
  Core::MathUtils::RotationMatrix3x3 beam_axis_rotation =
      MathUtils::rotation_matrix_from_direction(params.average_px, params.average_py, params.average_pz);

  size_t random_seed = std::stoull(IoUtils::get_required(config, "random_seed"));
  std::vector<Electron> beam;

  if (params.polar_lattice) {
    // Deterministic equal-area rings (IoUtils::polar_lattice_geometry); the generator only draws momenta. Sheets sit
    // at z_m = (m - (M-1)/2) L/M, moved by the modulation if it is on, each carrying sheet_weights' charge (1/M, or
    // the Gaussian profile's share when sigma_z > 0); M = 1 is the flat disk at z = 0 with unit weights.
    auto [rings, spacing] = IoUtils::polar_lattice_geometry(params.radius, params.lattice_spacing);
    const size_t sheets = params.stack_sheets;
    const std::vector<double> sheet_z = sheet_positions(params);
    const std::vector<double> sheet_w = sheet_weights(params);
    beam.reserve(sheets * IoUtils::polar_lattice_site_count(rings));
    std::mt19937 gen(random_seed);
    for (size_t m = 0; m < sheets; ++m) {
      const double z = sheet_z[m];
      const double charge_weight = sheet_w[m];
      beam.push_back(place_electron(params, beam_axis_rotation, 0.0, 0.0, z, tau_0, d_tau, N_tau, gen, charge_weight));
      for (int j = 1; j <= rings; ++j) {
        const double r = j * spacing;
        const long nodes = std::lround(2.0 * MathUtils::pi * j);
        for (long i = 0; i < nodes; ++i) {
          const double phi = 2.0 * MathUtils::pi * static_cast<double>(i) / static_cast<double>(nodes);
          beam.push_back(place_electron(params, beam_axis_rotation, r * std::cos(phi), r * std::sin(phi), z, tau_0,
                                        d_tau, N_tau, gen, charge_weight));
        }
      }
    }
    return beam;
  }

  size_t num_particles = std::stoull(IoUtils::get_required(config, "beam_particle_count"));
  beam.reserve(num_particles);  // Allocate memory upfront for speed

  // Fixed seed for reproducibility in physics testing; one generator shared
  // across every particle so the whole beam is a single deterministic
  // random sequence.
  std::mt19937 gen(random_seed);
  for (size_t p = 0; p < num_particles; ++p) {
    beam.push_back(generate_electron(params, beam_axis_rotation, tau_0, d_tau, N_tau, gen));
  }

  return beam;
}

}  // namespace Core::Particle
