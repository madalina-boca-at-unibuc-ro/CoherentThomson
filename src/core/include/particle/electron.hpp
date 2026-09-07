#pragma once
#include <iosfwd>
#include <vector>

#include "../laser/laser_field.hpp"
#include "../math_utils/math_utils.hpp"

namespace Core::Particle {

class Electron {
public:
  // One recorded point along the trajectory: the full state at a given
  // proper time
  struct State {
    Core::MathUtils::RealFourVector position;
    Core::MathUtils::RealFourVector momentum;
    double tau;
    // du^mu/dtau = (q_0/m_0) F^{mu nu}(x) u_nu, the electron's exact 4-acceleration at this state
    // (equivalently compute_derivative(position, momentum, laser).dmomentum / m_0). Filled in by
    // update_state()/compute_trajectory() -- see their comments below; zero-initialized (not the
    // true acceleration) for any state read back before compute_trajectory() has run.
    Core::MathUtils::RealFourVector acceleration;
  };

  // Constructs the electron from its initial 4-position and 4-velocity
  // (u^\mu). The current state is seeded from the same values; proper time
  // starts at initial_tau. The full trajectory (including this initial
  // state) is recorded on every subsequent set_state() call, with storage
  // reserved upfront for N_tau points.
  Electron(const Core::MathUtils::RealFourVector& initial_position,
           const Core::MathUtils::RealFourVector& initial_momentum, const double initial_tau, const double d_tau,
           const size_t N_tau);
  ~Electron();

  // Advances the current state by one proper-time step d_tau: reads the
  // current position/momentum, takes a 4th-order Runge-Kutta step of the
  // equations of motion
  //   dx^\mu/d\tau = p^\mu / m_0
  //   dp^\mu/d\tau = (q_0/m_0) F^{\mu\nu}(x) p_\nu
  // (with F the laser's Faraday tensor, and dp^\mu/d\tau the electron's
  // 4-acceleration up to the m_0 factor), overwrites the current
  // position/momentum with the result, and accumulates d_tau into the
  // proper time. The RK4 stage k1 (the exact derivative at the state being
  // stepped away from, i.e. trajectory.back() on entry) is reused, at no
  // extra cost, to fill in that state's acceleration -- see State::acceleration.
  void update_state(double d_tau, const Core::Laser::LaserField& laser);

  // Drives the electron through the rest of its N_tau-point trajectory: the
  // initial state is already recorded (at construction), so this calls
  // update_state(d_tau, laser) N_tau-1 more times to reach the full N_tau
  // recorded states, each with State::acceleration filled in via the k1 reuse
  // described above. The very last state has no following update_state() call
  // to supply that k1, so its acceleration is filled with one extra explicit
  // compute_derivative() call here (the only additional Faraday-tensor
  // evaluation this whole scheme costs).
  void compute_trajectory(const Core::Laser::LaserField& laser);

  // Initial state accessors
  const Core::MathUtils::RealFourVector& get_initial_position() const { return initial_position; }
  const Core::MathUtils::RealFourVector& get_initial_momentum() const { return initial_momentum; }

  // Current state accessors
  const Core::MathUtils::RealFourVector& get_position() const { return current_position; }
  const Core::MathUtils::RealFourVector& get_momentum() const { return current_momentum; }

  double get_proper_time() const { return proper_time; }
  double get_tau_0() const { return tau_0; }
  double get_d_tau() const { return d_tau; }
  size_t get_N_tau() const { return N_tau; }

  // The recorded trajectory
  const std::vector<State>& get_trajectory() const { return trajectory; }

  // Writes one row of "tau position(4) momentum(4) acceleration(4)" (no
  // trailing separator/newline) to stream. index == 0 always works (it reads
  // the initial position/momentum directly, not trajectory[0], though its
  // acceleration still comes from trajectory[0]); index > 0 requires
  // index < get_trajectory().size().
  void export_state(const size_t& index, std::ostream& stream) const;

private:
  double tau_0;  // initial value of the proper time
  double d_tau;  // tau_step in the numerical integration
  size_t N_tau;  // number of steps in the numerical integration
  // these are stored in each electron.

  // Phase-space derivative (dx^\mu/d\tau, dp^\mu/d\tau) at a given
  // position/momentum, per the Lorentz force law using the laser's Faraday
  // tensor at that position
  struct Derivative {
    Core::MathUtils::RealFourVector dposition;
    Core::MathUtils::RealFourVector dmomentum;
  };
  Derivative compute_derivative(const Core::MathUtils::RealFourVector& position,
                                const Core::MathUtils::RealFourVector& momentum,
                                const Core::Laser::LaserField& laser) const;

  // Initial state: 4-position x^\mu, 4-velocity u^\mu
  Core::MathUtils::RealFourVector initial_position;
  Core::MathUtils::RealFourVector initial_momentum;

  // Current state, mutated in-place by set_state() as the trajectory is
  // integrated
  Core::MathUtils::RealFourVector current_position;
  Core::MathUtils::RealFourVector current_momentum;

  // Proper time elapsed since the initial state; accumulates on every
  // set_state() call
  double proper_time = 0.0;

  // Full recorded trajectory, one State per set_state() call plus the
  // initial state
  std::vector<State> trajectory;
};

}  // namespace Core::Particle
