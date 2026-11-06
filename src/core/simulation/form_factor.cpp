#include "../include/simulation/form_factor.hpp"

#include <cmath>
#include <stdexcept>

#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Simulation {

FormFactor compute_form_factor(const std::vector<Particle::Electron>& electron_beam, const Laser::LaserField& laser,
                               const Detector::Detector_2D& detector) {
  if (electron_beam.empty()) {
    throw std::runtime_error("compute_form_factor: the electron beam is empty");
  }
  const double c = PhysUtils::AtomicUnits::c;
  const double mc = PhysUtils::AtomicUnits::m_0 * c;

  double W = 0.0;
  double px = 0.0, py = 0.0, pz = 0.0;
  for (const Particle::Electron& electron : electron_beam) {
    const MathUtils::RealFourVector& p = electron.get_initial_momentum();
    const double w = electron.get_charge_weight();
    W += w;
    px += w * p[1];
    py += w * p[2];
    pz += w * p[3];
  }
  px /= W;
  py /= W;
  pz /= W;
  MathUtils::RealFourVector P(std::sqrt(px * px + py * py + pz * pz + mc * mc), px, py, pz);

  // Same k1 and dressed momentum as Simulation::init_simulation_parameters, but n2 is the detector's actual normal
  // (including any tilt) and P is the generated beam's mean rather than the configured average_p.
  MathUtils::RealFourVector k_L = laser.get_unity_n() * laser.get_omega() / c;
  MathUtils::RealFourVector dressed = PhysUtils::dressed_momentum(P, k_L, laser.get_a0());
  const std::array<double, 3>& n = detector.get_normal_direction();
  MathUtils::RealFourVector n_s(1.0, n[0], n[1], n[2]);
  double omega_s_over_c = PhysUtils::non_linear_Thomson_formula(k_L, dressed, n_s, 1);
  MathUtils::RealFourVector k_s = n_s * omega_s_over_c;
  // The phase a translated electron's field picks up, measured against the solver (theory/dev_notes.md): shifting it
  // along its free worldline keeps the laser phase fixed, which brings in the undressed K_s.P / K_L.P.
  MathUtils::RealFourVector q = k_s - k_L * (MathUtils::contract(P, k_s) / MathUtils::contract(P, k_L));

  std::complex<double> A(0.0, 0.0);
  for (const Particle::Electron& electron : electron_beam) {
    double phase = MathUtils::dot3(q, electron.get_initial_position());
    A += electron.get_charge_weight() * std::complex<double>(std::cos(phase), std::sin(phase));
  }

  return {electron_beam.size(), W, P, dressed, omega_s_over_c * c, k_L, k_s, q, k_s - k_L, A / W};
}

}  // namespace Core::Simulation
