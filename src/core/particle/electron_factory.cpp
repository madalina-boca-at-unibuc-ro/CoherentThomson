#include "../include/particle/electron_factory.hpp"

#include <algorithm>
#include <iostream>
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

  // Convert cylindrical coordinates to Cartesian coordinates (cylinder height axis along canonical Oz).
  return place_electron(params, beam_axis_rotation, r * std::cos(phi), r * std::sin(phi), z, tau_0, d_tau, N_tau, gen);
}

ModulationGeometry modulation_geometry(const CylinderBeamParams& params, const Core::Laser::LaserField& laser,
                                       const std::array<double, 3>& detector_normal) {
  const double c = PhysUtils::AtomicUnits::c;
  const double mc = PhysUtils::AtomicUnits::m_0 * c;
  const double two_pi = 2.0 * MathUtils::pi;
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
  ModulationGeometry geometry;
  geometry.q_axial = q[1] * axis[0] + q[2] * axis[1] + q[3] * axis[2];
  const double q_squared = q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  geometry.q_transverse = std::sqrt(std::max(0.0, q_squared - geometry.q_axial * geometry.q_axial));
  if (std::abs(geometry.q_axial) < 1e-6 * std::sqrt(q_squared)) {
    throw std::invalid_argument("beam_modulation_zeta: q is perpendicular to the beam axis, so sheets along the axis "
                                "can't be matched to it");
  }

  const IoUtils::CylinderBeamParams::Modulation& modulation = params.modulation;
  const size_t sheets = params.stack_sheets;
  geometry.lambda_1 = two_pi / std::abs(geometry.q_axial);
  geometry.lambda_mod = modulation.period_ratio * geometry.lambda_1;
  geometry.K = two_pi / geometry.lambda_mod;
  geometry.spacing = geometry.lambda_mod / static_cast<double>(modulation.sheets_per_period);
  geometry.length = static_cast<double>(sheets) * geometry.spacing;
  geometry.spill = 0.0;
  geometry.sheet_z.reserve(sheets);
  for (size_t m = 0; m < sheets; ++m) {
    const double z_uniform = (static_cast<double>(m) - 0.5 * static_cast<double>(sheets - 1)) * geometry.spacing;
    const double z = z_uniform + modulation.zeta / geometry.K * std::sin(geometry.K * z_uniform + modulation.phase);
    geometry.sheet_z.push_back(z);
    geometry.spill = std::max(geometry.spill, std::abs(z) - geometry.length / 2.0);
  }
  return geometry;
}

std::vector<ModulationHarmonic> modulation_harmonics(const CylinderBeamParams& params,
                                                     const ModulationGeometry& geometry, size_t N_min,
                                                     size_t N_harmonics) {
  const IoUtils::CylinderBeamParams::Modulation& modulation = params.modulation;
  const double two_pi = 2.0 * MathUtils::pi;
  const double sign = geometry.q_axial > 0.0 ? 1.0 : -1.0;
  const bool whole_periods = params.stack_sheets % modulation.sheets_per_period == 0;
  std::vector<ModulationHarmonic> harmonics;
  for (size_t N = N_min; N < N_min + N_harmonics; ++N) {
    const double n = static_cast<double>(N) * modulation.period_ratio;
    std::complex<double> sum(0.0, 0.0);
    for (double z : geometry.sheet_z) {
      sum += std::polar(1.0, sign * n * geometry.K * z);
    }
    ModulationHarmonic harmonic{N, n, sum / static_cast<double>(geometry.sheet_z.size()), std::nullopt};
    // The midpoint rule on this periodic integrand converges fast once it has more points than n (1 + zeta), so it
    // stands in for (-1)^n J_n(n zeta) without needing a Bessel function from the standard library.
    if (whole_periods && n == std::round(n)) {
      const size_t points = 8 * static_cast<size_t>(std::ceil(n * (1.0 + modulation.zeta)) + 64.0);
      std::complex<double> cold(0.0, 0.0);
      for (size_t i = 0; i < points; ++i) {
        const double phi_0 = (static_cast<double>(i) + 0.5) * two_pi / static_cast<double>(points);
        const double phi = phi_0 + modulation.zeta * std::sin(phi_0);
        cold += std::polar(1.0, sign * n * (phi - modulation.phase));
      }
      harmonic.cold = cold / static_cast<double>(points);
    }
    harmonics.push_back(harmonic);
  }
  return harmonics;
}

std::vector<std::string> modulation_warnings(const CylinderBeamParams& params, size_t N_max) {
  const IoUtils::CylinderBeamParams::Modulation& modulation = params.modulation;
  const double S = static_cast<double>(modulation.sheets_per_period);
  const double M = static_cast<double>(params.stack_sheets);
  std::vector<std::string> warnings;
  if (params.stack_sheets % modulation.sheets_per_period != 0) {
    std::ostringstream text;
    text << "the stack holds " << M / S << " periods, not a whole number: b doesn't reduce to the cold-beam value";
    warnings.push_back(text.str());
  }
  const double recommended = 3.0 * static_cast<double>(N_max) * modulation.period_ratio;
  if (S < recommended) {
    std::ostringstream text;
    text << modulation.sheets_per_period
         << " sheets per period is below 3 N_max Lambda_mod / Lambda_1 = " << recommended
         << ": moved sheets can give a wrong b at the highest harmonics; see the log's sheet sum against the cold-beam "
            "value";
    warnings.push_back(text.str());
  }
  return warnings;
}

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config, const double tau_0, const double d_tau,
                                             const size_t N_tau, const Core::Laser::LaserField& laser,
                                             const std::array<double, 3>& detector_normal) {
  CylinderBeamParams params = IoUtils::parse_cylinder_beam_params(config);

  // Aligns the cylinder's height axis with the beam's own mean momentum direction (average_px/py/pz), so the
  // finite-height extent lies along the direction the beam actually moves in rather than always along canonical
  // Oz. Falls back to the identity if the mean momentum is (numerically) zero.
  Core::MathUtils::RotationMatrix3x3 beam_axis_rotation =
      MathUtils::rotation_matrix_from_direction(params.average_px, params.average_py, params.average_pz);

  size_t random_seed = std::stoull(IoUtils::get_required(config, "random_seed"));
  std::vector<Electron> beam;

  if (params.polar_lattice) {
    // Deterministic equal-area rings (IoUtils::polar_lattice_geometry); the generator only draws momenta. Sheets sit
    // at z_m = (m - (M-1)/2) d, each carrying 1/M of the charge; M = 1 is the flat disk at z = 0 with unit weights.
    // A bunched stack then moves them (modulation_geometry).
    auto [rings, spacing] = IoUtils::polar_lattice_geometry(params.radius, params.lattice_spacing);
    const size_t sheets = params.stack_sheets;
    const double charge_weight = 1.0 / static_cast<double>(sheets);
    std::vector<double> sheet_z;
    if (params.modulation.active) {
      sheet_z = modulation_geometry(params, laser, detector_normal).sheet_z;
      const size_t N_max = IoUtils::get_number_of_harmonics_min(config) + IoUtils::get_number_of_harmonics(config) - 1;
      for (const std::string& warning : modulation_warnings(params, N_max)) {
        std::cerr << "Warning: " << warning << "\n";
      }
    }
    beam.reserve(sheets * IoUtils::polar_lattice_site_count(rings));
    std::mt19937 gen(random_seed);
    for (size_t m = 0; m < sheets; ++m) {
      const double z = params.modulation.active
                           ? sheet_z[m]
                           : (static_cast<double>(m) - 0.5 * static_cast<double>(sheets - 1)) * params.stack_spacing;
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
