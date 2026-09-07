#include "../include/particle/electron.hpp"

#include <iostream>
#include <ostream>

#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Particle {

Electron::Electron(const Core::MathUtils::RealFourVector& initial_position,
                   const Core::MathUtils::RealFourVector& initial_momentum, const double initial_tau,
                   const double d_tau, const size_t N_tau)
    : tau_0(initial_tau),
      d_tau(d_tau),
      N_tau(N_tau),
      initial_position(initial_position),
      initial_momentum(initial_momentum),
      current_position(initial_position),
      current_momentum(initial_momentum),
      proper_time(initial_tau) {
  trajectory.reserve(N_tau);
  trajectory.push_back({initial_position, initial_momentum, initial_tau, Core::MathUtils::RealFourVector{}});
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

  // k1 = compute_derivative(x0, p0, laser) is the *exact* phase-space derivative at (x0, p0), i.e.
  // the state trajectory.back() already holds on entry (either the initial state, or the state
  // pushed by the previous update_state() call) -- so k1.dmomentum / m_0 is exactly that state's
  // 4-acceleration (State::acceleration), obtained here for free rather than via an extra
  // Faraday-tensor evaluation.
  trajectory.back().acceleration = k1.dmomentum / Core::PhysUtils::AtomicUnits::m_0;

  current_position = x0 + (d_tau / 6.0) * (k1.dposition + 2.0 * k2.dposition + 2.0 * k3.dposition + k4.dposition);
  current_momentum = p0 + (d_tau / 6.0) * (k1.dmomentum + 2.0 * k2.dmomentum + 2.0 * k3.dmomentum + k4.dmomentum);

  proper_time += d_tau;

  trajectory.push_back({current_position, current_momentum, proper_time, Core::MathUtils::RealFourVector{}});
}

void Electron::compute_trajectory(const Core::Laser::LaserField& laser) {
  for (size_t i = 0; i < N_tau - 1; ++i) {
    update_state(d_tau, laser);
  }
  // The loop above fills every state's acceleration from the *next* state's RK4 k1 (see
  // update_state); the final state has no following call to supply that k1, so it's filled in here
  // with one extra, otherwise-unneeded derivative evaluation.
  trajectory.back().acceleration =
      compute_derivative(current_position, current_momentum, laser).dmomentum / Core::PhysUtils::AtomicUnits::m_0;
}

void Electron::export_state(const size_t& index, std::ostream& stream) const {
  const MathUtils::RealFourVector& position = index == 0 ? initial_position : trajectory[index].position;
  const MathUtils::RealFourVector& momentum = index == 0 ? initial_momentum : trajectory[index].momentum;
  double tau = index == 0 ? tau_0 : trajectory[index].tau;
  const MathUtils::RealFourVector& acceleration = trajectory[index].acceleration;

  stream << tau << " " << position[0] << " " << position[1] << " " << position[2] << " " << position[3] << " ";
  stream << momentum[0] << " " << momentum[1] << " " << momentum[2] << " " << momentum[3] << " ";
  stream << acceleration[0] << " " << acceleration[1] << " " << acceleration[2] << " " << acceleration[3];
}

}  // namespace Core::Particle
