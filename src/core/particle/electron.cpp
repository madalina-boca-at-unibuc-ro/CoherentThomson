#include "../include/particle/electron.hpp"

#include <iostream>
#include <ostream>
#include <stdexcept>

#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Particle {

Electron::Electron(const Core::MathUtils::RealFourVector& initial_position,
                   const Core::MathUtils::RealFourVector& initial_momentum, const double initial_tau,
                   const double d_tau, const size_t N_tau, const bool store_trajectory_in)
    : tau_0(initial_tau),
      d_tau(d_tau),
      N_tau(N_tau),
      initial_position(initial_position),
      initial_momentum(initial_momentum),
      current_position(initial_position),
      current_momentum(initial_momentum),
      proper_time(initial_tau),
      store_trajectory(store_trajectory_in) {
  if (store_trajectory) {
    trajectory.reserve(N_tau);
    trajectory.push_back({initial_position, initial_momentum, initial_tau});
  }
}

Electron::~Electron() = default;

Electron::Derivative Electron::compute_derivative(const Core::MathUtils::RealFourVector& position,
                                                  const Core::MathUtils::RealFourVector& momentum,
                                                  const Core::Laser::LaserField& laser) const {
  // dx^\mu/d\tau = p^\mu / m_0
  Core::MathUtils::RealFourVector dposition = momentum / Core::PhysUtils::AtomicUnits::m_0;

  // dp^\mu/d\tau = (q_0/m_0) F^{\mu\nu}(x) p_\nu
  Core::Laser::FaradayTensor F = laser.get_faraday_tensor(position);
  Core::MathUtils::RealFourVector dmomentum =
      (Core::PhysUtils::AtomicUnits::q_0 / Core::PhysUtils::AtomicUnits::m_0) * Core::MathUtils::contract(F, momentum);

  return {dposition, dmomentum};
}

void Electron::update_state(double d_tau, const Core::Laser::LaserField& laser) {
  const Core::MathUtils::RealFourVector& x0 = current_position;
  const Core::MathUtils::RealFourVector& p0 = current_momentum;

  // Classic 4th-order Runge-Kutta step over the coupled (position, momentum) state
  Derivative k1 = compute_derivative(x0, p0, laser);
  Derivative k2 = compute_derivative(x0 + (d_tau / 2.0) * k1.dposition, p0 + (d_tau / 2.0) * k1.dmomentum, laser);
  Derivative k3 = compute_derivative(x0 + (d_tau / 2.0) * k2.dposition, p0 + (d_tau / 2.0) * k2.dmomentum, laser);
  Derivative k4 = compute_derivative(x0 + d_tau * k3.dposition, p0 + d_tau * k3.dmomentum, laser);

  current_position = x0 + (d_tau / 6.0) * (k1.dposition + 2.0 * k2.dposition + 2.0 * k3.dposition + k4.dposition);
  current_momentum = p0 + (d_tau / 6.0) * (k1.dmomentum + 2.0 * k2.dmomentum + 2.0 * k3.dmomentum + k4.dmomentum);

  proper_time += d_tau;

  if (store_trajectory) {
    trajectory.push_back({current_position, current_momentum, proper_time});
  }
}

void Electron::compute_trajectory(const Core::Laser::LaserField& laser) {
  if (!store_trajectory) {
    throw std::runtime_error(
        "Electron::compute_trajectory: called on an electron constructed with store_trajectory_in = false");
  }
  for (size_t i = 0; i < N_tau - 1; ++i) {
    update_state(d_tau, laser);
  }
}

void Electron::export_state(const size_t& index, std::ostream& stream) const {
  const MathUtils::RealFourVector& position = index == 0 ? initial_position : trajectory[index].position;
  const MathUtils::RealFourVector& momentum = index == 0 ? initial_momentum : trajectory[index].momentum;
  double tau = index == 0 ? tau_0 : trajectory[index].tau;

  stream << tau << " " << position[0] << " " << position[1] << " " << position[2] << " " << position[3] << " ";
  stream << momentum[0] << " " << momentum[1] << " " << momentum[2] << " " << momentum[3];
}

}  // namespace Core::Particle
