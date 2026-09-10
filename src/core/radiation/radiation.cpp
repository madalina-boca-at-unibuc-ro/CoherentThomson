#include "../include/radiation/radiation.hpp"

#include <algorithm>
#include <cmath>

namespace Core::Radiation {

// ComplexBivector aliases std::array, whose only associated namespace (for ADL) is std -- pull
// MathUtils::operator+= into scope explicitly so `bivector += bivector` below resolves.
using Core::MathUtils::operator+=;

namespace {

// One element of the antisymmetric bivector n_R0^\mu u^\nu - n_R0^\nu u^\mu, for alpha < beta.
// F[beta][alpha] is always -F[alpha][beta] and the diagonal is identically zero, so only
// the 6 independent upper-triangle elements of the (otherwise 4x4) tensor are ever
// evaluated below -- never the full 16.
inline double bivector_element(const MathUtils::RealFourVector& n, const MathUtils::RealFourVector& u, size_t alpha,
                               size_t beta) {
  return n[alpha] * u[beta] - n[beta] * u[alpha];
}

// Accumulates amp_long*term_long / amp_short*term_short into slot `index` (0..5, in the fixed
// (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) order -- see MathUtils::ComplexBivector) of
// `long_bivector`/`short_bivector`. term_long/term_short are frequency-independent geometric
// factors precomputed once per (tau, screen point) and shared across every frequency -- the
// simplified form (see below) shares one bare n^alpha u^beta - n^beta u^alpha factor between both,
// but the direct form's long-range term needs a different combined (n, u, w) factor, hence the two
// separate parameters rather than one shared `term`. Only the 6 independent upper-triangle
// elements of the corresponding Faraday tensor are ever accumulated -- the antisymmetric lower
// triangle and zero diagonal are filled in once, after all contributions (from every tau,
// electron, and thread) have been summed, via MathUtils::unpack_bivector, called once in
// Simulation::run_simulation on the final reduced field.
inline void add_bivector_term(MathUtils::ComplexBivector& long_bivector, MathUtils::ComplexBivector& short_bivector,
                              size_t index, double term_long, double term_short, const MathUtils::Complex& amp_long,
                              const MathUtils::Complex& amp_short) {
  long_bivector[index] += amp_long * term_long;
  short_bivector[index] += amp_short * term_short;
}

// Accumulates one signed endpoint contribution into slot `index` of `boundary_bivector` -- see
// "Form 2's boundary term F_b" below. Kept as its own one-line helper (rather than folded into
// add_bivector_term above) since it's only ever called from the two endpoint branches
// (i_tau == 0 / i_tau == N_tau - 1), not from every tau like add_bivector_term.
inline void add_boundary_term(MathUtils::ComplexBivector& boundary_bivector, size_t index, double term,
                              const MathUtils::Complex& amp_boundary) {
  boundary_bivector[index] += amp_boundary * term;
}

// ---- Integrand construction: PREFACT * exp(i * freq * phase_argument) ----
//
// The per-(tau, screen point, frequency) integrand factors into an exponential phase term
// (radiation_phase_argument) and two PREFACT terms (long_range_prefactor/short_range_prefactor)
// that each multiply the shared geometric bivector term. The current PREFACT formulas were
// obtained by integrating the standard radiation integral by parts, and are the piece most
// likely to change if a different derivation is adopted (see CLAUDE.md's "OPEN VALIDATION GAP"
// note) -- isolated into their own functions below for exactly that reason, so a new formula can
// be dropped in here without touching the phase/loop-structure code around them.

// The "rest of the exponent" in PREFACT * exp(i * N*omega * phase_argument): frequency itself is
// applied by the caller, not here, so the fast evenly-spaced-frequency phase-factor recurrence in
// compute_radiation's i_freq loop stays valid regardless of how this function or the PREFACT
// terms below change. Shared verbatim by both the simplified and direct forms (see
// theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Common notation" section) -- the two
// forms differ only in their PREFACT terms, never in this phase.
inline double radiation_phase_argument(const MathUtils::RealFourVector& x, double R) { return x[0] + R; }

// Long-range PREFACT term (multiplies the shared bivector term n_R0^alpha u^beta - n_R0^beta
// u^alpha). Frequency-dependent: proportional to freq/R. `inv_R` (rather than R) is taken as a
// parameter so callers that sweep many frequencies at fixed (tau, screen point) can precompute
// the division once.
inline MathUtils::Complex long_range_prefactor(double freq, double inv_R) {
  return MathUtils::Complex{0.0, -freq * inv_R};
}

// Short-range PREFACT term (multiplies the same shared bivector term). Frequency-independent in
// the current formula. Negative sign: d/dtau(1/|R_0|) = +(n_R0.u)/|R_0|^2 (three-vector dot), so
// the integration-by-parts term -integral(f' * B) that produces this piece carries an overall
// minus -- see theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Sign bug (fixed)" note
// under Form 2.
inline double short_range_prefactor(const MathUtils::RealFourVector& n_R0, const MathUtils::RealFourVector& u, double R,
                                    double n_R0_contract_u) {
  return -MathUtils::dot3(n_R0, u) / (R * R * n_R0_contract_u);
}

// ---- Boundary term (theory/FT_Faraday_tensor-direct_and_simplified_forms.md, "Form 2's boundary
// term F_b") -- only meaningful for the simplified formula, and only nonzero at the two ends of the
// trajectory's finite tau range. ----
//
// The integration by parts that produces F_l/F_s from Jackson's compact form also produces a
// boundary term F_b = [e^{ik*Phi(tau)}/|R_0(tau)| * T^{alpha beta}(tau)]_{tau_min}^{tau_max}, where
// T^{alpha beta} = (n_R0^alpha u^beta - n_R0^beta u^alpha)/(n_R0.u) is the same tensor factor F_l/F_s
// share. This is only negligible for an integral over all of tau in (-infinity, infinity); since
// compute_radiation integrates each electron over a *finite* recorded trajectory, F_b is generally
// significant and is not optional -- see the theory doc's "On-axis F^03 cancellation" section, which
// shows it can be as large as F_l itself near the beam axis. F_b's prefactor shares the shared
// bivector term T^{alpha beta} = term_{alpha beta}/n_R0_contract_u with F_s's short-range term
// above, so only n_R0_contract_u's reciprocal is new here; the endpoint sign (+1 at tau_max, -1 at
// tau_min) is applied by the caller, which is the only place that knows a given tau is an endpoint.
inline double boundary_prefactor(double inv_R, double n_R0_contract_u) { return inv_R / n_R0_contract_u; }

// ---- Direct form (theory/FT_Faraday_tensor-direct_and_simplified_forms.md, Form 1) ----
//
// The direct FT of the closed-form Lienard-Wiechert field, as opposed to the simplified form's
// integration-by-parts derivation above: neither PREFACT term below carries an explicit frequency
// factor (the only omega-dependence is in the shared exp(i*omega*phase_argument) factor, applied
// by the caller), and the long-range term's geometric factor is a genuine combination of the
// (n_R0, u) and (n_R0, w) bivectors, not the bare (n_R0, u) bivector alone -- see
// direct_long_range_tensor_term below. Selected via the "radiation_formula"="direct" config key.

// Direct-form long-range PREFACT term: 1 / (R * (n_R0.u)^2). Denominator power is 2, not 3 -- an
// earlier version of this formula (and of the theory doc it implements) dropped the Jacobian
// dt/d(tau) = (u.n_R0)/c when changing the outer integration variable from t to tau, which left an
// extra spurious power of (u.n_R0) in the denominator; see the theory doc's Form 1 "Jacobian bug
// (fixed)" note.
inline double long_range_prefactor_direct(double R, double n_R0_contract_u) {
  return 1.0 / (R * n_R0_contract_u * n_R0_contract_u);
}

// Direct-form short-range PREFACT term: 1 / (R^2 * (n_R0.u)^2). No longer carries an explicit c^2
// factor: the same Jacobian fix that drops the denominator power to 2 (see
// long_range_prefactor_direct above) also folds what used to be this term's separate c^2 into the
// overall e/(4 pi eps0 c^2) normalization constant now shared uniformly by F_l and F_s -- applied
// once, for every formula, via Simulation::run_simulation's general_factor, rather than baked in
// here per formula as before.
inline double short_range_prefactor_direct(double R, double n_R0_contract_u) {
  return 1.0 / (R * R * n_R0_contract_u * n_R0_contract_u);
}

// Direct-form long-range geometric factor: (n_R0.u)*(n_R0^alpha w^beta - n_R0^beta w^alpha) -
// (n_R0.w)*(n_R0^alpha u^beta - n_R0^beta u^alpha), replacing the simplified form's bare (n_R0, u)
// bivector for the long-range term only -- the short-range term still uses the bare (n_R0, u)
// bivector (bivector_element(n_R0, u, alpha, beta)) in both forms.
inline double direct_long_range_tensor_term(const MathUtils::RealFourVector& n_R0, const MathUtils::RealFourVector& u,
                                            const MathUtils::RealFourVector& w, double n_R0_contract_u,
                                            double n_R0_contract_w, size_t alpha, size_t beta) {
  return n_R0_contract_u * bivector_element(n_R0, w, alpha, beta) -
         n_R0_contract_w * bivector_element(n_R0, u, alpha, beta);
}

}  // namespace

void compute_radiation(Particle::Electron& electron, const Laser::LaserField& laser,
                       const std::vector<double>& frequencies_list, const Detector::Detector_2D& detector,
                       Simulation::PackedRadiationField& field, bool use_direct_formula) {
  electron.compute_trajectory(laser);

  size_t N_tau = electron.get_N_tau();
  size_t N_d = detector.get_total_points();
  size_t N_freq = frequencies_list.size();
  const std::vector<Particle::Electron::State>& trajectory = electron.get_trajectory();
  double d_tau = electron.get_d_tau();

  // Simulation::init_simulation_parameters builds frequencies_list as an arithmetic progression in
  // both its modes: N_freq consecutive harmonics N_harmonics_min, N_harmonics_min+1, ... of a
  // single fundamental (non_linear_Thomson_formula is exactly linear in the harmonic index, so
  // frequencies_list[i] == frequencies_list[0] + i*step for a constant step -- see simulation.cpp),
  // or (dense_frequency_spectrum=true) a plain linear scan between omega_min/omega_max. Either way,
  // exp(i*phase*frequencies_list[i]) is just exp(i*phase*frequencies_list[0]) times
  // exp(i*phase*step) raised to the i-th power, so every entry's phase factor can be obtained from
  // two cos/sin evaluations plus cheap complex multiplications instead of N_freq separate
  // transcendental evaluations -- checked once here (not per trajectory/screen point) so a future
  // change to non-evenly-spaced frequencies safely falls back to the direct per-frequency
  // evaluation below rather than silently computing the wrong phase. Note this is unrelated to
  // which harmonic index the list starts at (N_harmonics_min): the recurrence only needs constant
  // spacing, not that frequencies_list[0] itself be the fundamental.
  double step = N_freq > 1 ? frequencies_list[1] - frequencies_list[0] : 0.0;
  bool frequencies_are_evenly_spaced = N_freq > 0;
  for (size_t i = 2; frequencies_are_evenly_spaced && i < N_freq; ++i) {
    double expected = frequencies_list[0] + static_cast<double>(i) * step;
    frequencies_are_evenly_spaced = std::abs(frequencies_list[i] - expected) <= 1e-9 * std::abs(expected);
  }

  // Per-frequency accumulators for one screen point at a time, reused (and zeroed) across
  // screen points rather than reallocated every iteration. Packed as ComplexBivector (6 complex
  // elements) rather than the full 4x4 ComplexFourTensor, since only the 6 independent
  // upper-triangle Faraday elements are ever accumulated here -- keeps this array small enough to
  // stay cache/register-resident across the N_tau sweep below even for a large N_freq.
  std::vector<MathUtils::ComplexBivector> local_long(N_freq);
  std::vector<MathUtils::ComplexBivector> local_short(N_freq);
  // Boundary-term accumulator (see boundary_prefactor's doc comment above) -- always allocated but
  // only ever written to for the simplified formula, and only at the two trajectory endpoints
  // (i_tau == 0 / i_tau == N_tau - 1), so it stays zero for the direct formula and for every
  // interior tau either way.
  std::vector<MathUtils::ComplexBivector> local_boundary(N_freq);

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
    std::fill(local_boundary.begin(), local_boundary.end(), MathUtils::ComplexBivector{});

    for (size_t i_tau = 0; i_tau < N_tau; i_tau++) {
      const MathUtils::RealFourVector& x = trajectory[i_tau].position;
      const MathUtils::RealFourVector& u = trajectory[i_tau].momentum;

      // R_0 = x_0 - r_0(tau) (observation point minus particle position), matching
      // theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Common notation" section --
      // n_R0 = R_0/|R_0| points from the emitting particle to the observer, not the reverse.
      MathUtils::RealFourVector diff = detector_point - x;

      // n_R0 is built in place from diff, which also yields R (the spatial norm) without
      // recomputing the sqrt a second time.
      MathUtils::RealFourVector n_R0 = diff;
      double R = MathUtils::create_unit_light_like_vector_in_place(n_R0);

      double n_R0_contract_u = MathUtils::contract(n_R0, u);

      // The bare (n_R0, u) bivector terms depend only on n_R0 and u, not on frequency, so they're
      // computed once per (tau, screen point) here rather than once per frequency below. Used
      // directly as the short-range term in both formulas, and as the long-range term in the
      // simplified formula (see long_term01..23 below for the direct formula's own long-range
      // factor).
      double term01 = bivector_element(n_R0, u, 0, 1);
      double term02 = bivector_element(n_R0, u, 0, 2);
      double term03 = bivector_element(n_R0, u, 0, 3);
      double term12 = bivector_element(n_R0, u, 1, 2);
      double term13 = bivector_element(n_R0, u, 1, 3);
      double term23 = bivector_element(n_R0, u, 2, 3);

      // amp_short_0/inv_R (simplified formula) or prefactor_l_direct/prefactor_s_direct (direct
      // formula) are the tau/screen-point-level pieces of the PREFACT terms that don't depend on
      // frequency, so they're computed once per tau here rather than once per frequency below.
      // long_term01..23 is the long-range geometric factor: the bare (n_R0, u) bivector for the
      // simplified formula, or direct_long_range_tensor_term's (n_R0, u, w) combination for the
      // direct formula -- see radiation.hpp's use_direct_formula doc comment and
      // theory/FT_Faraday_tensor-direct_and_simplified_forms.md.
      double amp_short_0 = 0.0;
      double inv_R = 1.0 / R;
      double prefactor_l_direct = 0.0;
      double prefactor_s_direct = 0.0;
      double long_term01 = term01, long_term02 = term02, long_term03 = term03;
      double long_term12 = term12, long_term13 = term13, long_term23 = term23;

      if (use_direct_formula) {
        const MathUtils::RealFourVector& w = trajectory[i_tau].acceleration;
        double n_R0_contract_w = MathUtils::contract(n_R0, w);
        prefactor_l_direct = long_range_prefactor_direct(R, n_R0_contract_u);
        prefactor_s_direct = short_range_prefactor_direct(R, n_R0_contract_u);
        long_term01 = direct_long_range_tensor_term(n_R0, u, w, n_R0_contract_u, n_R0_contract_w, 0, 1);
        long_term02 = direct_long_range_tensor_term(n_R0, u, w, n_R0_contract_u, n_R0_contract_w, 0, 2);
        long_term03 = direct_long_range_tensor_term(n_R0, u, w, n_R0_contract_u, n_R0_contract_w, 0, 3);
        long_term12 = direct_long_range_tensor_term(n_R0, u, w, n_R0_contract_u, n_R0_contract_w, 1, 2);
        long_term13 = direct_long_range_tensor_term(n_R0, u, w, n_R0_contract_u, n_R0_contract_w, 1, 3);
        long_term23 = direct_long_range_tensor_term(n_R0, u, w, n_R0_contract_u, n_R0_contract_w, 2, 3);
      } else {
        amp_short_0 = short_range_prefactor(n_R0, u, R, n_R0_contract_u);
      }

      // Boundary term only exists for the simplified formula (see boundary_prefactor's doc comment
      // above), and only contributes at the two ends of the finite tau range: -1 at tau_min
      // (i_tau == 0), +1 at tau_max (i_tau == N_tau - 1). Both conditions hold simultaneously when
      // N_tau == 1, in which case the two signed contributions are equal and opposite and cancel
      // exactly -- the correct limit for a zero-width integration interval -- so no special-casing
      // is needed beyond evaluating both independently. amp_boundary_sign is folded in here (rather
      // than at the accumulation site) so the i_freq loop below only ever needs to check "is this
      // an endpoint at all", not which one.
      bool is_lower_endpoint = !use_direct_formula && i_tau == 0;
      bool is_upper_endpoint = !use_direct_formula && i_tau == N_tau - 1;
      double boundary_weight = 0.0;
      if (is_lower_endpoint) boundary_weight -= boundary_prefactor(inv_R, n_R0_contract_u);
      if (is_upper_endpoint) boundary_weight += boundary_prefactor(inv_R, n_R0_contract_u);
      bool has_boundary_contribution = is_lower_endpoint || is_upper_endpoint;

      // Trapezoidal quadrature weight approximating the continuous tau-integral that F_l/F_s
      // are defined as (theory/FT_Faraday_tensor-direct_and_simplified_forms.md's F_l/F_s formulas
      // are literally \int d\tau ... over the electron's finite recorded trajectory). Electron
      // trajectories are recorded at a fixed RK4 step (Electron::get_d_tau), so this reduces to the
      // standard uniform-grid trapezoidal rule: full d_tau weight at every interior tau, half weight
      // at the two endpoints. Deliberately NOT applied to the boundary term below: F_b is an exact
      // antiderivative evaluation at the two endpoints, not itself a tau-sum, so it needs no
      // quadrature weight (verified against the independent Python reference implementation, which
      // treats it the same way). Previously missing entirely, which made every long-range/
      // short-range contribution scale with N_tau (number of trajectory samples) instead of
      // converging to the physical integral as the trajectory is resolved more finely.
      double tau_weight = d_tau;
      if (N_tau > 1 && (i_tau == 0 || i_tau == N_tau - 1)) tau_weight *= 0.5;

      double phase_base = radiation_phase_argument(x, R);

      // cexp, frequencies_list[0]'s phase factor, and cexp_step, the constant spacing's phase
      // factor, are each computed via one cos/sin evaluation (std::polar(1, theta) ==
      // exp(i*theta) without the generic complex-exp path's extra real-exponential evaluation).
      // When frequencies_are_evenly_spaced holds, every subsequent entry's phase factor is
      // obtained by multiplying by cexp_step again rather than by another transcendental
      // evaluation. Guarded by N_freq > 0 since frequencies_list[0] would otherwise be an
      // out-of-bounds read. Both formulas share the exact same phase (see
      // radiation_phase_argument's doc comment), so this construction is unaffected by
      // use_direct_formula.
      MathUtils::Complex cexp = N_freq > 0 ? std::polar(1.0, phase_base * frequencies_list[0]) : MathUtils::Complex{};
      MathUtils::Complex cexp_step = N_freq > 1 ? std::polar(1.0, phase_base * step) : MathUtils::Complex{};

      for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
        if (!frequencies_are_evenly_spaced) {
          cexp = std::polar(1.0, phase_base * frequencies_list[i_freq]);
        }
        // The direct formula's PREFACT terms carry no explicit frequency factor (only the shared
        // phase cexp above does), unlike the simplified formula's long-range term -- see
        // long_range_prefactor_direct/short_range_prefactor_direct's doc comments.
        MathUtils::Complex amp_long = use_direct_formula ? prefactor_l_direct * cexp
                                                         : long_range_prefactor(frequencies_list[i_freq], inv_R) * cexp;
        MathUtils::Complex amp_short = use_direct_formula ? prefactor_s_direct * cexp : amp_short_0 * cexp;
        amp_long *= tau_weight;
        amp_short *= tau_weight;

        add_bivector_term(local_long[i_freq], local_short[i_freq], 0, long_term01, term01, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 1, long_term02, term02, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 2, long_term03, term03, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 3, long_term12, term12, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 4, long_term13, term13, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 5, long_term23, term23, amp_long, amp_short);

        // Only the two endpoint tau values ever reach this (has_boundary_contribution is false,
        // hence amp_boundary needs never be formed, for every interior tau) -- negligible added
        // cost relative to the N_tau-wide long/short accumulation above.
        if (has_boundary_contribution) {
          MathUtils::Complex amp_boundary = boundary_weight * cexp;
          add_boundary_term(local_boundary[i_freq], 0, term01, amp_boundary);
          add_boundary_term(local_boundary[i_freq], 1, term02, amp_boundary);
          add_boundary_term(local_boundary[i_freq], 2, term03, amp_boundary);
          add_boundary_term(local_boundary[i_freq], 3, term12, amp_boundary);
          add_boundary_term(local_boundary[i_freq], 4, term13, amp_boundary);
          add_boundary_term(local_boundary[i_freq], 5, term23, amp_boundary);
        }

        if (frequencies_are_evenly_spaced) cexp *= cexp_step;
      }
    }

    for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
      field.field[i_freq][i_d].long_range += local_long[i_freq];
      field.field[i_freq][i_d].short_range += local_short[i_freq];
      field.field[i_freq][i_d].boundary += local_boundary[i_freq];
    }
  }
}

}  // namespace Core::Radiation
