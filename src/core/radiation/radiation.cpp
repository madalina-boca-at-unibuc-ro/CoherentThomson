#include "../include/radiation/radiation.hpp"

#include <algorithm>
#include <cmath>

namespace Core::Radiation {

// ComplexBivector aliases std::array, whose only associated namespace (for ADL) is std -- pull
// MathUtils::operator+= into scope explicitly so `bivector += bivector` below resolves.
using Core::MathUtils::operator+=;

namespace {

// One element of the antisymmetric bivector n^\mu u^\nu - n^\nu u^\mu, for alpha < beta.
// F[beta][alpha] is always -F[alpha][beta] and the diagonal is identically zero, so only
// the 6 independent upper-triangle elements of the (otherwise 4x4) tensor are ever
// evaluated below -- never the full 16.
inline double bivector_element(const MathUtils::RealFourVector& n, const MathUtils::RealFourVector& u, size_t alpha,
                               size_t beta) {
  return n[alpha] * u[beta] - n[beta] * u[alpha];
}

// Accumulates amp_long/amp_short * term into slot `index` (0..5, in the fixed
// (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) order -- see MathUtils::ComplexBivector) of
// `long_bivector`/`short_bivector`, where `term` is the frequency-independent geometric factor
// n^alpha u^beta - n^beta u^alpha, precomputed once per (tau, screen point) and shared across
// every frequency. Only the 6 independent upper-triangle elements of the corresponding Faraday
// tensor are ever accumulated -- the antisymmetric lower triangle and zero diagonal are filled in
// once, after all contributions (from every tau, electron, and thread) have been summed, via
// MathUtils::unpack_bivector, called once in Simulation::run_simulation on the final reduced
// field.
inline void add_bivector_term(MathUtils::ComplexBivector& long_bivector, MathUtils::ComplexBivector& short_bivector,
                              size_t index, double term, const MathUtils::Complex& amp_long,
                              const MathUtils::Complex& amp_short) {
  long_bivector[index] += amp_long * term;
  short_bivector[index] += amp_short * term;
}

}  // namespace

void compute_radiation(Particle::Electron& electron, const Laser::LaserField& laser,
                       const std::vector<double>& frequencies_list, const Detector::Detector_2D& detector,
                       Simulation::PackedRadiationField& field) {
  electron.compute_trajectory(laser);

  size_t N_tau = electron.get_N_tau();
  size_t N_d = detector.get_total_points();
  size_t N_freq = frequencies_list.size();
  const std::vector<Particle::Electron::State>& trajectory = electron.get_trajectory();

  // Simulation::init_simulation_parameters builds frequencies_list as the first N_freq harmonics
  // of a single fundamental (frequencies_list[i] == (i+1) * frequencies_list[0] exactly, since
  // only the harmonic index N varies across non_linear_Thomson_formula calls -- see
  // simulation.cpp). When that holds, exp(i*phase*frequencies_list[i]) is just
  // exp(i*phase*frequencies_list[0]) raised to the (i+1)-th power, so every harmonic's phase
  // factor can be obtained from one cos/sin evaluation plus cheap complex multiplications instead
  // of N_freq separate transcendental evaluations -- checked once here (not per trajectory/screen
  // point) so a future change to non-harmonic frequencies (e.g. honoring omega_min/omega_max,
  // see CLAUDE.md) safely falls back to the direct per-frequency evaluation below rather than
  // silently computing the wrong phase.
  bool frequencies_are_harmonics = N_freq > 0;
  for (size_t i = 1; frequencies_are_harmonics && i < N_freq; ++i) {
    double expected = static_cast<double>(i + 1) * frequencies_list[0];
    frequencies_are_harmonics = std::abs(frequencies_list[i] - expected) <= 1e-9 * std::abs(expected);
  }

  // Per-frequency accumulators for one screen point at a time, reused (and zeroed) across
  // screen points rather than reallocated every iteration. Packed as ComplexBivector (6 complex
  // elements) rather than the full 4x4 ComplexFourTensor, since only the 6 independent
  // upper-triangle Faraday elements are ever accumulated here -- keeps this array small enough to
  // stay cache/register-resident across the N_tau sweep below even for a large N_freq.
  std::vector<MathUtils::ComplexBivector> local_long(N_freq);
  std::vector<MathUtils::ComplexBivector> local_short(N_freq);

  // Here frequencies are actually divided by c; see PhysUtils.hpp
  //
  // Screen points (i_d) are the outer loop and trajectory points (i_tau) the inner loop,
  // deliberately the reverse of the natural "for each tau, splat onto every screen point"
  // order: `field` is the full [N_freq][N_d] output grid (tens of MB for a fine detector),
  // while `trajectory` is comparatively tiny (N_tau states, already fully materialized above).
  // Every trajectory point contributes to every (screen point, frequency) pair regardless of
  // loop order, so the total arithmetic is identical either way -- but with i_tau innermost,
  // each screen point's full tau-sum lands in a couple of small local tensors that stay hot in
  // cache/registers, and the large `field` array is only touched once per (i_d, i_freq) pair
  // after the tau loop completes, instead of once per (i_tau, i_d, i_freq) triple. That avoids
  // re-sweeping the entire multi-MB `field` array on every one of the N_tau trajectory steps,
  // which otherwise dominates runtime as pure memory traffic once N_tau*N_d*N_freq is large.
  for (size_t i_d = 0; i_d < N_d; i_d++) {
    const MathUtils::RealFourVector detector_point = detector.get_point(i_d);
    std::fill(local_long.begin(), local_long.end(), MathUtils::ComplexBivector{});
    std::fill(local_short.begin(), local_short.end(), MathUtils::ComplexBivector{});

    for (size_t i_tau = 0; i_tau < N_tau; i_tau++) {
      const MathUtils::RealFourVector& x = trajectory[i_tau].position;
      const MathUtils::RealFourVector& u = trajectory[i_tau].momentum;

      MathUtils::RealFourVector diff = x - detector_point;

      // n0 is built in place from diff, which also yields R (the spatial norm) without
      // recomputing the sqrt a second time.
      MathUtils::RealFourVector n0 = diff;
      double R = MathUtils::create_unit_light_like_vector_in_place(n0);

      double n0_contract_u = MathUtils::contract(n0, u) * R;
      double n0_dot3_u = MathUtils::dot3(n0, u);
      // amp_short_0 (real) and amp_long_0 = -i/R (purely imaginary) are kept as plain doubles
      // rather than Complex(amp_short_0, 0) / Complex(0, -1/R): below, amp_short = amp_short_0 *
      // cexp and amp_long = amp_long_0 * freq * cexp then reduce to a handful of real multiplies
      // directly from cexp's own real/imaginary parts, instead of two generic
      // complex-times-complex multiplications per frequency.
      double amp_short_0 = n0_dot3_u / (R * R * n0_contract_u);
      double inv_R = 1.0 / R;

      // The geometric bivector terms depend only on n0 and u, not on frequency, so they're
      // computed once per (tau, screen point) here rather than once per frequency below.
      double term01 = bivector_element(n0, u, 0, 1);
      double term02 = bivector_element(n0, u, 0, 2);
      double term03 = bivector_element(n0, u, 0, 3);
      double term12 = bivector_element(n0, u, 1, 2);
      double term13 = bivector_element(n0, u, 1, 3);
      double term23 = bivector_element(n0, u, 2, 3);

      double phase_base = x[0] + R;

      // cexp_1, the fundamental's phase factor, is computed via one cos/sin evaluation
      // (std::polar(1, theta) == exp(i*theta) without the generic complex-exp path's extra
      // real-exponential evaluation). When frequencies_are_harmonics holds, every subsequent
      // harmonic's phase factor is obtained by multiplying by cexp_1 again rather than by another
      // transcendental evaluation. Guarded by N_freq > 0 since frequencies_list[0] would
      // otherwise be an out-of-bounds read.
      MathUtils::Complex cexp_1 = N_freq > 0 ? std::polar(1.0, phase_base * frequencies_list[0]) : MathUtils::Complex{};
      MathUtils::Complex cexp = cexp_1;

      for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
        if (!frequencies_are_harmonics) {
          cexp = std::polar(1.0, phase_base * frequencies_list[i_freq]);
        }
        // amp_long = (-i * freq/R) * cexp and amp_short = amp_short_0 * cexp, expanded directly
        // from cexp = c + i*s (see the amp_short_0/inv_R comment above).
        double c = cexp.real();
        double s = cexp.imag();
        double w = frequencies_list[i_freq] * inv_R;
        MathUtils::Complex amp_long{w * s, -w * c};
        MathUtils::Complex amp_short{amp_short_0 * c, amp_short_0 * s};

        add_bivector_term(local_long[i_freq], local_short[i_freq], 0, term01, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 1, term02, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 2, term03, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 3, term12, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 4, term13, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 5, term23, amp_long, amp_short);

        if (frequencies_are_harmonics) cexp *= cexp_1;
      }
    }

    for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
      field.field[i_freq][i_d].long_range += local_long[i_freq];
      field.field[i_freq][i_d].short_range += local_short[i_freq];
    }
  }
}

}  // namespace Core::Radiation
