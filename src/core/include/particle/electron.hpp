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
  };

  // Constructs the electron from its initial 4-position and 4-velocity
  // (u^\mu). The current state is seeded from the same values; proper time
  // starts at initial_tau. When store_trajectory_in is true, the full
  // trajectory (including this initial state) is recorded on every
  // subsequent set_state() call, with storage reserved upfront for
  // reserve_steps points; when false (the default), no trajectory storage is
  // allocated at all.
  Electron(const Core::MathUtils::RealFourVector& initial_position,
           const Core::MathUtils::RealFourVector& initial_momentum, const double initial_tau, const double d_tau,
           const size_t N_tau, const bool store_trajectory_in = false);
  ~Electron();

  // Advances the current state by one proper-time step d_tau: reads the
  // current position/momentum, takes a 4th-order Runge-Kutta step of the
  // equations of motion
  //   dx^\mu/d\tau = p^\mu / m_0
  //   dp^\mu/d\tau = (q_0/m_0) F^{\mu\nu}(x) p_\nu
  // (with F the laser's Faraday tensor, and dp^\mu/d\tau the electron's
  // 4-acceleration up to the m_0 factor), overwrites the current
  // position/momentum with the result, and accumulates d_tau into the
  // proper time.
  void update_state(double d_tau, const Core::Laser::LaserField& laser);

  // Drives the electron through the rest of its N_tau-point trajectory: the
  // initial state is already recorded (at construction, when
  // store_trajectory_in was true), so this calls update_state(d_tau, laser)
  // N_tau-1 more times to reach the full N_tau recorded states. Throws
  // std::runtime_error if the electron was constructed with
  // store_trajectory_in = false, since there is no trajectory to complete.
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
  bool get_store_trajectory() const { return store_trajectory; }

  // The recorded trajectory (empty unless store_trajectory_in was true at
  // construction)
  const std::vector<State>& get_trajectory() const { return trajectory; }

  // Writes one row of "tau position(4) momentum(4)" (no trailing
  // separator/newline) to stream. index == 0 always works (it reads the
  // initial state directly, not trajectory[0]), even when no trajectory was
  // recorded; index > 0 requires get_store_trajectory() and
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

  // Whether set_state() appends to trajectory; false by default so large
  // runs pay no storage cost
  bool store_trajectory = false;

  // Full recorded trajectory, one State per set_state() call plus the
  // initial state; only populated when store_trajectory is true
  std::vector<State> trajectory;
};

}  // namespace Core::Particle
