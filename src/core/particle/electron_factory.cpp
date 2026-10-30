#include "../include/particle/electron_factory.hpp"

namespace Core::Particle {

namespace {

// Places one electron at (x, y, z) in the cylinder's canonical frame and draws its momentum: translation by
// beam_center, rotation onto the mean momentum, then Gaussian momenta about average_p (advances gen).
Electron place_electron(const CylinderBeamParams& params, const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation,
                        double x, double y, double z, double tau_0, double d_tau, size_t N_tau, std::mt19937& gen) {
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

  return Electron(position, momentum, tau_0, d_tau, N_tau);
}

}  // namespace

Electron generate_electron(const CylinderBeamParams& params, const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation,
                           double tau_0, double d_tau, size_t N_tau, std::mt19937& gen) {
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
  return place_electron(params, beam_axis_rotation, r * std::cos(phi), r * std::sin(phi), z, tau_0, d_tau, N_tau,
                        gen);
}

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config, const double tau_0, const double d_tau,
                                             const size_t N_tau) {
  CylinderBeamParams params = IoUtils::parse_cylinder_beam_params(config);

  // Aligns the cylinder's height axis with the beam's own mean momentum direction (average_px/py/pz), so the
  // finite-height extent lies along the direction the beam actually moves in rather than always along canonical
  // Oz. Falls back to the identity if the mean momentum is (numerically) zero.
  Core::MathUtils::RotationMatrix3x3 beam_axis_rotation =
      MathUtils::rotation_matrix_from_direction(params.average_px, params.average_py, params.average_pz);

  size_t random_seed = std::stoull(IoUtils::get_required(config, "random_seed"));
  std::vector<Electron> beam;

  if (params.polar_lattice) {
    // Deterministic equal-area rings (IoUtils::polar_lattice_geometry); the generator only draws momenta.
    auto [rings, spacing] = IoUtils::polar_lattice_geometry(params.radius, params.lattice_spacing);
    std::mt19937 gen(random_seed);
    beam.push_back(place_electron(params, beam_axis_rotation, 0.0, 0.0, 0.0, tau_0, d_tau, N_tau, gen));
    for (int j = 1; j <= rings; ++j) {
      const double r = j * spacing;
      const long nodes = std::lround(2.0 * MathUtils::pi * j);
      for (long i = 0; i < nodes; ++i) {
        const double phi = 2.0 * MathUtils::pi * static_cast<double>(i) / static_cast<double>(nodes);
        beam.push_back(place_electron(params, beam_axis_rotation, r * std::cos(phi), r * std::sin(phi), 0.0, tau_0,
                                      d_tau, N_tau, gen));
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
