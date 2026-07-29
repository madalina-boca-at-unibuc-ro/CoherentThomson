#include "../include/particle/electron_factory.hpp"

namespace Core::Particle {

Electron generate_electron(const CylinderBeamParams& params, const Core::MathUtils::RealFourTensor& rotation_matrix,
                           const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation, double tau_0, double d_tau,
                           size_t N_tau, std::mt19937& gen, bool store_trajectory) {
  std::uniform_real_distribution<double> dis_uniform(0.0, 1.0);
  std::uniform_real_distribution<double> dis_z(-params.height / 2.0, params.height / 2.0);

  // Area-preserving radial distribution transformation
  double r = params.radius * std::sqrt(dis_uniform(gen));
  double phi = dis_uniform(gen) * 2.0 * MathUtils::pi;
  double z = dis_z(gen);

  // Convert cylindrical coordinates to Cartesian coordinates (cylinder height axis along canonical Oz), then
  // apply the beam's global translation before either rotation below, so the offset itself is expressed in the
  // same canonical frame as beam_cylinder_radius/height and gets carried along with the beam through both
  // rotations (rather than being a fixed lab-frame shift applied afterwards).
  double x = r * std::cos(phi) + params.center_x;
  double y = r * std::sin(phi) + params.center_y;
  z += params.center_z;

  // Reorient the cylinder so its height axis matches the beam's own mean momentum direction, before the whole
  // beam is rotated together with the laser below.
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

  // Rotate the initial electron distribution and momenta together with the laser's actual direction.
  position = MathUtils::contract(rotation_matrix, position);
  momentum = MathUtils::contract(rotation_matrix, momentum);

  return Electron(position, momentum, tau_0, d_tau, N_tau, store_trajectory);
}

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config,
                                             const Core::MathUtils::RealFourTensor rotation_matrix, const double tau_0,
                                             const double d_tau, const size_t N_tau) {
  CylinderBeamParams params = IoUtils::parse_cylinder_beam_params(config);

  // Aligns the cylinder's height axis with the beam's own mean momentum direction (average_px/py/pz), so the
  // finite-height extent lies along the direction the beam actually moves in rather than always along canonical
  // Oz. Falls back to the identity if the mean momentum is (numerically) zero.
  Core::MathUtils::RotationMatrix3x3 beam_axis_rotation =
      MathUtils::rotation_matrix_from_direction(params.average_px, params.average_py, params.average_pz);

  size_t num_particles = std::stoull(IoUtils::get_required(config, "beam_particle_count"));
  size_t random_seed = std::stoull(IoUtils::get_required(config, "random_seed"));

  // Trajectory storage: off by default (zero overhead for large runs); when
  // enabled, each electron reserves exactly the number of RK4 steps the
  // laser's own phase grid implies (laser_NT points per cycle times the
  // pulse length in cycles), so no reallocation happens during the run.
  bool store_trajectory = IoUtils::get_required(config, "store_trajectory") == "true";

  std::vector<Electron> beam;
  beam.reserve(num_particles);  // Allocate memory upfront for speed

  // Fixed seed for reproducibility in physics testing; one generator shared
  // across every particle so the whole beam is a single deterministic
  // random sequence.
  std::mt19937 gen(random_seed);
  for (size_t p = 0; p < num_particles; ++p) {
    beam.push_back(
        generate_electron(params, rotation_matrix, beam_axis_rotation, tau_0, d_tau, N_tau, gen, store_trajectory));
  }

  return beam;
}

}  // namespace Core::Particle
