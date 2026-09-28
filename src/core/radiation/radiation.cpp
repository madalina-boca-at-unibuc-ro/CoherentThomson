#include "../include/radiation/radiation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include "../include/phys_utils/phys_utils.hpp"

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
// factors precomputed once per (tau, screen point) and shared across every frequency. They differ
// in both forms: the simplified form's short-range term uses the (s, u) bivector with s = (0, n_R0)
// (short_range_bivector_element), and the direct form's long-range term a combined (n, u, w)
// factor (direct_long_range_tensor_term), hence two separate parameters rather than one `term`. Only the 6 independent
// upper-triangle elements of the corresponding Faraday tensor are ever accumulated -- the antisymmetric lower triangle
// and zero diagonal are filled in once, after all contributions (from every tau, electron, and thread) have been
// summed, via MathUtils::unpack_bivector, called once in Simulation::run_simulation on the final reduced field.
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
// that each multiply their own geometric bivector term. The current PREFACT formulas were
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

// Short-range PREFACT term: 1/R^2, frequency-independent (theory/
// FT_Faraday_tensor-direct_and_simplified_forms-v2.md, Form 2, Eq. IV.1.2.15). Multiplies the
// short-range bivector short_range_bivector_element(n_R0, u, ...) below, NOT the (n_R0, u) bivector
// F_l/F_b use. This replaces the v1 kernel -(n_R0.u)_3 / (R^2 (n_R0.u)) * (n_R0, u)-bivector, which came
// from treating Jackson's fixed-event d/dtau as a derivative along the retarded observation path
// and so lost finite-distance terms. The v2 kernel instead follows from differentiating the
// Fourier-transformed Lienard-Wiechert potential A^beta ~ int dtau e^{ik Phi} u^beta/|R_0|: the
// spatial gradient of the 1/|R_0| amplitude is the only source of 1/|R_0|^2, hence s = (0, n_R0).
inline double short_range_prefactor(double inv_R) { return inv_R * inv_R; }

// One element of the simplified form's short-range bivector s^alpha u^beta - s^beta u^alpha, for
// alpha < beta, with s = (0, n_R0_x, n_R0_y, n_R0_z) (theory doc v2, Form 2). Equal to
// bivector_element(n_R0, u, ...) for spatial (alpha, beta); for (0, beta) it is -n_R0^beta u^0
// instead of u^beta - n_R0^beta u^0, since s^0 = 0. beta is never 0 because alpha < beta.
inline double short_range_bivector_element(const MathUtils::RealFourVector& n, const MathUtils::RealFourVector& u,
                                           size_t alpha, size_t beta) {
  double s_alpha = alpha == 0 ? 0.0 : n[alpha];
  return s_alpha * u[beta] - n[beta] * u[alpha];
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
// shows it can be as large as F_l itself near the beam axis. F_b multiplies the same (n_R0, u)
// bivector as F_l (T^{alpha beta} = term_{alpha beta}/n_R0_contract_u); unchanged by the v2
// correction of the theory doc, which only replaced F_s. The endpoint sign (+1 at tau_max, -1 at
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

// Direct-form short-range PREFACT term: c^2 / (R^2 * (n_R0.u)^2). The theory doc's F_s carries a
// bare q/(2 pi) (atomic units), while Simulation::run_simulation's general_factor applies
// q/(2 pi c^2) uniformly to every term, so the c^2 must be restored here. Commit be8cccd dropped it
// (reading the doc's "F_s loses its separate c^2" as removing it from the code as well), which made
// the direct F_s ~c^2 = 1.9e4 times too small -- found by the simplified-vs-direct near-field
// cross-check that validated the v2 simplified F_s.
inline double short_range_prefactor_direct(double R, double n_R0_contract_u) {
  return PhysUtils::AtomicUnits::c * PhysUtils::AtomicUnits::c / (R * R * n_R0_contract_u * n_R0_contract_u);
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

namespace {

// Shared by both exact and long-distance paths: frequencies_list is normally an arithmetic
// progression (see compute_radiation_exact's comment on the recurrence), so each tau's phase factors
// can be built by recurrence from two std::polar calls instead of one per frequency.
bool are_evenly_spaced(const std::vector<double>& frequencies_list, double step) {
  size_t N_freq = frequencies_list.size();
  bool evenly_spaced = N_freq > 0;
  for (size_t i = 2; evenly_spaced && i < N_freq; ++i) {
    double expected = frequencies_list[0] + static_cast<double>(i) * step;
    evenly_spaced = std::abs(frequencies_list[i] - expected) <= 1e-9 * std::abs(expected);
  }
  return evenly_spaced;
}

// Exact formulas ("simplified"/"direct"): theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md.
void compute_radiation_exact(const Particle::Electron& electron, const std::vector<double>& frequencies_list,
                             const Detector::Detector_2D& detector, Simulation::PackedRadiationField& field,
                             bool use_direct_formula) {
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
  bool frequencies_are_evenly_spaced = are_evenly_spaced(frequencies_list, step);

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
      // as the short-range term in the direct formula, and as the long-range and boundary terms in
      // the simplified formula (see long_term01..23 below for the direct formula's own long-range
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
      // Short-range geometric factor: the bare (n_R0, u) bivector for the direct formula, the
      // (s, u) bivector (s = (0, n_R0)) for the simplified formula -- see
      // short_range_bivector_element.
      double short_term01 = term01, short_term02 = term02, short_term03 = term03;
      double short_term12 = term12, short_term13 = term13, short_term23 = term23;

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
        amp_short_0 = short_range_prefactor(inv_R);
        short_term01 = short_range_bivector_element(n_R0, u, 0, 1);
        short_term02 = short_range_bivector_element(n_R0, u, 0, 2);
        short_term03 = short_range_bivector_element(n_R0, u, 0, 3);
        short_term12 = short_range_bivector_element(n_R0, u, 1, 2);
        short_term13 = short_range_bivector_element(n_R0, u, 1, 3);
        short_term23 = short_range_bivector_element(n_R0, u, 2, 3);
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

        add_bivector_term(local_long[i_freq], local_short[i_freq], 0, long_term01, short_term01, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 1, long_term02, short_term02, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 2, long_term03, short_term03, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 3, long_term12, short_term12, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 4, long_term13, short_term13, amp_long, amp_short);
        add_bivector_term(local_long[i_freq], local_short[i_freq], 5, long_term23, short_term23, amp_long, amp_short);

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

// Long-distance formulas ("long_distance_simplified"/"long_distance_direct"):
// theory/long_distance_direct_simplified_coding_guide.md, sections 1-6, via the two-projection method.
//
// Per (electron, screen point) the geometry is frozen at the electron's first sample r(tau_m):
// x_0 = x - r(tau_m), n_0 = x_0/|x_0|, and a fixed transverse basis e_1, e_2 (e_1 x e_2 = n_0). Per tau
// only two complex scalars per frequency are accumulated (p_1, p_2), and the tensor is rebuilt once per
// frequency as p_a (e_a^alpha n_0^beta - n_0^alpha e_a^beta) / |x_0|, general_factor being applied by
// run_simulation as for the exact formulas. The full phase k(|x_0| + r^0 - n_0.r_0) is split into the
// constant k(|x_0| + r^0(tau_m)), applied once per frequency at reconstruction, and the varying
// k(r^0 - r^0(tau_m) - n_0.r_0) inside the tau loop. This is only for precision: |x_0| ~ 1e9 a.u. would
// otherwise dominate the argument and cost ~1e-9 rad of rounding per tau, amplified by the tau-sum
// cancellation at weak harmonics. The in-loop phase stays frequency-independent, so the evenly-spaced
// recurrence applies. Only the spatial position is shifted by r(tau_m); the time is not reset (its constant
// part is in the reconstruction factor), so relative phases between electrons are kept.
//
// GCC's SLP vectorizer miscompiles this function under -O3 -march=native (seen with GCC 15.2): the
// complex accumulation in the i_freq loop below then returns wrong sums at every frequency after the
// first (h2 off by ~100x, h3 by ~1e4x), while -O0/-O2/-O3 without -march=native, or -O3 -march=native
// with -fno-tree-slp-vectorize or -ffp-contract=off, all agree. UBSan/ASan report nothing, and the
// per-tau inputs are correct, so this is taken as a compiler bug. SLP is therefore disabled for this
// function only (the exact path is unaffected). Re-check before removing: run long_distance_simplified
// with N_harmonics >= 2 and compare against a -O0 build.
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-tree-slp-vectorize")))
#endif
void compute_radiation_long_distance(const Particle::Electron& electron, const std::vector<double>& frequencies_list,
                                     const Detector::Detector_2D& detector, Simulation::PackedRadiationField& field,
                                     bool use_direct_formula) {
  size_t N_tau = electron.get_N_tau();
  size_t N_d = detector.get_total_points();
  size_t N_freq = frequencies_list.size();
  const std::vector<Particle::Electron::State>& trajectory = electron.get_trajectory();
  double d_tau = electron.get_d_tau();
  if (N_tau == 0 || N_freq == 0) return;

  double step = N_freq > 1 ? frequencies_list[1] - frequencies_list[0] : 0.0;
  bool frequencies_are_evenly_spaced = are_evenly_spaced(frequencies_list, step);

  const MathUtils::RealFourVector& reference = trajectory[0].position;  // r(tau_m); only [1..3] are used

  // bulk[i_freq][a]: the tau-integral part of p_a (simplified: without its ik factor, applied after the
  // tau loop; direct: the whole p_a). endpoint[i_freq][a]: the simplified endpoint bracket (with its sign).
  using Projections = std::array<MathUtils::Complex, 2>;
  std::vector<Projections> bulk(N_freq);
  std::vector<Projections> endpoint(N_freq);

  for (size_t i_d = 0; i_d < N_d; i_d++) {
    const MathUtils::RealFourVector detector_point = detector.get_point(i_d);
    std::fill(bulk.begin(), bulk.end(), Projections{});
    std::fill(endpoint.begin(), endpoint.end(), Projections{});

    std::array<double, 3> n0{detector_point[1] - reference[1], detector_point[2] - reference[2],
                             detector_point[3] - reference[3]};
    double x0_norm = std::sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
    if (!(x0_norm > 0.0)) {
      throw std::runtime_error("Radiation::compute_radiation: long-distance formula needs |x - r(tau_m)| > 0, "
                               "but screen point " +
                               std::to_string(i_d) + " coincides with an electron's start");
    }
    for (double& component : n0)
      component /= x0_norm;

    // Transverse basis (guide section 5): the Cartesian axis least aligned with n_0, projected
    // orthogonal to n_0 and normalized, then e_2 = n_0 x e_1.
    size_t axis = 0;
    for (size_t i = 1; i < 3; ++i) {
      if (std::abs(n0[i]) < std::abs(n0[axis])) axis = i;
    }
    std::array<double, 3> e1{};
    e1[axis] = 1.0;
    for (size_t i = 0; i < 3; ++i)
      e1[i] -= n0[axis] * n0[i];
    double e1_norm = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (double& component : e1)
      component /= e1_norm;
    std::array<double, 3> e2{n0[1] * e1[2] - n0[2] * e1[1], n0[2] * e1[0] - n0[0] * e1[2],
                             n0[0] * e1[1] - n0[1] * e1[0]};
    const std::array<std::array<double, 3>, 2> basis{e1, e2};

    for (size_t i_tau = 0; i_tau < N_tau; i_tau++) {
      const MathUtils::RealFourVector& x = trajectory[i_tau].position;
      const MathUtils::RealFourVector& u = trajectory[i_tau].momentum;

      double n0_dot_r0 = n0[0] * (x[1] - reference[1]) + n0[1] * (x[2] - reference[2]) + n0[2] * (x[3] - reference[3]);
      double n0_contract_u = u[0] - (n0[0] * u[1] + n0[1] * u[2] + n0[2] * u[3]);
      std::array<double, 2> e_dot_u{};
      for (size_t a = 0; a < 2; ++a)
        e_dot_u[a] = basis[a][0] * u[1] + basis[a][1] * u[2] + basis[a][2] * u[3];

      double tau_weight = d_tau;
      if (N_tau > 1 && (i_tau == 0 || i_tau == N_tau - 1)) tau_weight *= 0.5;

      // Frequency-independent geometric factor multiplying the phase in the tau-integral of p_a.
      // Simplified: e_a.u (times ik later). Direct: -[(n_0.u)(e_a.w) - (n_0.w)(e_a.u)]/(n_0.u)^2.
      std::array<double, 2> bulk_term{};
      if (use_direct_formula) {
        const MathUtils::RealFourVector& w = trajectory[i_tau].acceleration;
        double n0_contract_w = w[0] - (n0[0] * w[1] + n0[1] * w[2] + n0[2] * w[3]);
        for (size_t a = 0; a < 2; ++a) {
          double e_dot_w = basis[a][0] * w[1] + basis[a][1] * w[2] + basis[a][2] * w[3];
          bulk_term[a] =
              -tau_weight * (n0_contract_u * e_dot_w - n0_contract_w * e_dot_u[a]) / (n0_contract_u * n0_contract_u);
        }
      } else {
        for (size_t a = 0; a < 2; ++a)
          bulk_term[a] = tau_weight * e_dot_u[a];
      }

      // Simplified endpoint bracket -[(e_a.u)/(n_0.u) e^{ik Phi}]_{tau_m}^{tau_M}: +1 at tau_m, -1 at
      // tau_M after the leading minus; both apply (and cancel) when N_tau == 1. No quadrature weight.
      double endpoint_sign = 0.0;
      if (!use_direct_formula) {
        if (i_tau == 0) endpoint_sign += 1.0;
        if (i_tau == N_tau - 1) endpoint_sign -= 1.0;
      }
      std::array<double, 2> endpoint_term{};
      for (size_t a = 0; a < 2; ++a)
        endpoint_term[a] = endpoint_sign * e_dot_u[a] / n0_contract_u;
      bool has_endpoint_contribution = endpoint_sign != 0.0;

      double phase_base = (x[0] - reference[0]) - n0_dot_r0;
      MathUtils::Complex cexp = std::polar(1.0, phase_base * frequencies_list[0]);
      MathUtils::Complex cexp_step = N_freq > 1 ? std::polar(1.0, phase_base * step) : MathUtils::Complex{};

      for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
        if (!frequencies_are_evenly_spaced) cexp = std::polar(1.0, phase_base * frequencies_list[i_freq]);
        bulk[i_freq][0] += bulk_term[0] * cexp;
        bulk[i_freq][1] += bulk_term[1] * cexp;
        if (has_endpoint_contribution) {
          endpoint[i_freq][0] += endpoint_term[0] * cexp;
          endpoint[i_freq][1] += endpoint_term[1] * cexp;
        }
        if (frequencies_are_evenly_spaced) cexp *= cexp_step;
      }
    }

    // Reconstruction (guide section 6): sum_a p_a (e_a^alpha n_0^beta - n_0^alpha e_a^beta) / |x_0|, in
    // the packed (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) order, with e_a^0 = 0 and n_0^0 = 1.
    std::array<std::array<double, 6>, 2> basis_bivector{};
    for (size_t a = 0; a < 2; ++a) {
      const std::array<double, 3>& e = basis[a];
      basis_bivector[a] = {
          -e[0], -e[1], -e[2], e[0] * n0[1] - n0[0] * e[1], e[0] * n0[2] - n0[0] * e[2], e[1] * n0[2] - n0[1] * e[2]};
    }
    double inv_x0_norm = 1.0 / x0_norm;
    double constant_phase = x0_norm + reference[0];  // |x_0| + r^0(tau_m), taken out of the tau loop
    for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
      // e^{ik(|x_0| + r^0(tau_m))}/|x_0|; the simplified bulk p_a also carries ik (guide 5.2), while the
      // direct p_a is complete as accumulated (guide 5.1).
      MathUtils::Complex amplitude = std::polar(inv_x0_norm, frequencies_list[i_freq] * constant_phase);
      MathUtils::Complex bulk_factor =
          use_direct_formula ? amplitude : MathUtils::Complex{0.0, frequencies_list[i_freq]} * amplitude;
      Simulation::PackedFaraday& target = field.field[i_freq][i_d];
      for (size_t a = 0; a < 2; ++a) {
        MathUtils::Complex p_bulk = bulk_factor * bulk[i_freq][a];
        MathUtils::Complex p_endpoint = amplitude * endpoint[i_freq][a];
        for (size_t index = 0; index < 6; ++index) {
          target.long_range[index] += p_bulk * basis_bivector[a][index];
          target.boundary[index] += p_endpoint * basis_bivector[a][index];
        }
      }
    }
  }
}

// e^{i angle} for a small angle (|angle| <= kMaxSteppedAngle), by Taylor series: cos to angle^10 and sin
// to angle^11, so the truncation error is below ~1e-16 at the largest allowed angle.
inline void small_angle_phasor(double angle, double& cos_value, double& sin_value) {
  double a2 = angle * angle;
  cos_value =
      1.0 - a2 * (1.0 / 2) *
                (1.0 - a2 * (1.0 / 12) * (1.0 - a2 * (1.0 / 30) * (1.0 - a2 * (1.0 / 56) * (1.0 - a2 * (1.0 / 90)))));
  sin_value =
      angle * (1.0 - a2 * (1.0 / 6) *
                         (1.0 - a2 * (1.0 / 20) *
                                    (1.0 - a2 * (1.0 / 42) * (1.0 - a2 * (1.0 / 72) * (1.0 - a2 * (1.0 / 110))))));
}

constexpr double kMaxSteppedAngle = 0.25;  // polynomial range; larger angles are halved first
constexpr size_t kMaxHalvings = 6;         // up to |angle| = 16 rad per step; beyond that std::polar is used
constexpr size_t kAnchorInterval = 64;     // exact std::polar re-anchoring period, in tau steps
constexpr size_t kScreenBlock = 8;         // screen points stepped together (independent chains)

// Multiplies each phasor (re[l], im[l]) of a block by e^{i scale*increment[l]*2^halvings}: the polynomial
// at scale*increment[l] (|.| <= kMaxSteppedAngle), squared `halvings` times.
inline void advance_phasors(double* re, double* im, const double* increment, double scale, size_t halvings) {
  double c[kScreenBlock], s[kScreenBlock];
  for (size_t l = 0; l < kScreenBlock; ++l)
    small_angle_phasor(scale * increment[l], c[l], s[l]);
  for (size_t h = 0; h < halvings; ++h) {
    for (size_t l = 0; l < kScreenBlock; ++l) {
      double c_squared = c[l] * c[l] - s[l] * s[l];
      s[l] = 2.0 * c[l] * s[l];
      c[l] = c_squared;
    }
  }
  for (size_t l = 0; l < kScreenBlock; ++l) {
    double new_re = re[l] * c[l] - im[l] * s[l];
    im[l] = re[l] * s[l] + im[l] * c[l];
    re[l] = new_re;
  }
}

// Same formulas and output as compute_radiation_long_distance ("long_distance_*_approx" modes), but the
// phase factor is advanced along tau instead of evaluated from scratch. With n_0 fixed, the long-distance
// phase is linear in the trajectory, so between consecutive samples it changes by
// Delta Phi_j = Delta r^0_j - n_0.Delta r_j (differences of stored samples, hence accurate), and
// e^{ik Phi_{j+1}} = e^{ik Phi_j} e^{ik Delta Phi_j}, with the second factor from small_angle_phasor. Every
// kAnchorInterval steps (and if |k Delta Phi| exceeds even the halving range) the phasor is reset from the exact
// in-loop phase with std::polar, which bounds the accumulated rounding error. As in the unstepped function,
// the constant k(|x_0| + r^0(tau_m)) is applied at reconstruction, not in the loop. kScreenBlock screen points are
// processed together so their independent phasor chains overlap in the CPU (a single chain is
// latency-bound and barely beats std::polar). Frequencies use the same evenly-spaced recurrence as the other
// paths, with both of its phasors (base and step) stepped; a non-evenly-spaced list falls back to
// compute_radiation_long_distance.
//
// SLP is disabled here for the same reason as in compute_radiation_long_distance (GCC 15.2 miscompile).
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-tree-slp-vectorize")))
#endif
void compute_radiation_long_distance_stepped(const Particle::Electron& electron,
                                             const std::vector<double>& frequencies_list,
                                             const Detector::Detector_2D& detector,
                                             Simulation::PackedRadiationField& field, bool use_direct_formula) {
  size_t N_tau = electron.get_N_tau();
  size_t N_d = detector.get_total_points();
  size_t N_freq = frequencies_list.size();
  if (N_tau == 0 || N_freq == 0 || N_d == 0) return;

  double step = N_freq > 1 ? frequencies_list[1] - frequencies_list[0] : 0.0;
  if (!are_evenly_spaced(frequencies_list, step)) {
    compute_radiation_long_distance(electron, frequencies_list, detector, field, use_direct_formula);
    return;
  }
  const double k0 = frequencies_list[0];
  const bool multi_frequency = N_freq > 1;
  const std::vector<Particle::Electron::State>& trajectory = electron.get_trajectory();
  const double d_tau = electron.get_d_tau();
  const MathUtils::RealFourVector& reference = trajectory[0].position;  // r(tau_m)
  constexpr size_t L = kScreenBlock;

  // Screen-point-independent per-tau data, in structure-of-arrays form: r_0 = r - r(tau_m) and the
  // increments Delta r_j = r_j - r_{j-1} (index 0 unused), taken directly from the stored samples.
  std::vector<double> t0(N_tau), rx(N_tau), ry(N_tau), rz(N_tau), dt0(N_tau), drx(N_tau), dry(N_tau), drz(N_tau);
  for (size_t j = 0; j < N_tau; ++j) {
    const MathUtils::RealFourVector& x = trajectory[j].position;
    t0[j] = x[0] - reference[0];
    rx[j] = x[1] - reference[1];
    ry[j] = x[2] - reference[2];
    rz[j] = x[3] - reference[3];
    if (j > 0) {
      const MathUtils::RealFourVector& x_prev = trajectory[j - 1].position;
      dt0[j] = x[0] - x_prev[0];
      drx[j] = x[1] - x_prev[1];
      dry[j] = x[2] - x_prev[2];
      drz[j] = x[3] - x_prev[3];
    }
  }

  // Accumulators for one block, indexed [(i_freq * 2 + a) * L + l]: the tau-integral part of p_a (without
  // the simplified form's ik) and the simplified endpoint bracket, as in compute_radiation_long_distance.
  std::vector<double> bulk_re(N_freq * 2 * L), bulk_im(N_freq * 2 * L);
  std::vector<double> end_re(N_freq * 2 * L), end_im(N_freq * 2 * L);

  for (size_t block_begin = 0; block_begin < N_d; block_begin += L) {
    size_t block_size = std::min(L, N_d - block_begin);

    // Per-screen-point geometry; slots past block_size repeat the last point and are never stored.
    double D[L], nx[L], ny[L], nz[L], e1x[L], e1y[L], e1z[L], e2x[L], e2y[L], e2z[L];
    for (size_t l = 0; l < L; ++l) {
      size_t i_d = block_begin + std::min(l, block_size - 1);
      const MathUtils::RealFourVector detector_point = detector.get_point(i_d);
      std::array<double, 3> n0{detector_point[1] - reference[1], detector_point[2] - reference[2],
                               detector_point[3] - reference[3]};
      double x0_norm = std::sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
      if (!(x0_norm > 0.0)) {
        throw std::runtime_error("Radiation::compute_radiation: long-distance formula needs |x - r(tau_m)| > 0, "
                                 "but screen point " +
                                 std::to_string(i_d) + " coincides with an electron's start");
      }
      for (double& component : n0)
        component /= x0_norm;
      size_t axis = 0;
      for (size_t i = 1; i < 3; ++i) {
        if (std::abs(n0[i]) < std::abs(n0[axis])) axis = i;
      }
      std::array<double, 3> e1{};
      e1[axis] = 1.0;
      for (size_t i = 0; i < 3; ++i)
        e1[i] -= n0[axis] * n0[i];
      double e1_norm = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
      for (double& component : e1)
        component /= e1_norm;
      D[l] = x0_norm;
      nx[l] = n0[0];
      ny[l] = n0[1];
      nz[l] = n0[2];
      e1x[l] = e1[0];
      e1y[l] = e1[1];
      e1z[l] = e1[2];
      e2x[l] = n0[1] * e1[2] - n0[2] * e1[1];
      e2y[l] = n0[2] * e1[0] - n0[0] * e1[2];
      e2z[l] = n0[0] * e1[1] - n0[1] * e1[0];
    }

    std::fill(bulk_re.begin(), bulk_re.end(), 0.0);
    std::fill(bulk_im.begin(), bulk_im.end(), 0.0);
    std::fill(end_re.begin(), end_re.end(), 0.0);
    std::fill(end_im.begin(), end_im.end(), 0.0);

    // Base phasor e^{i k0 Phi} and frequency-step phasor e^{i step Phi}, carried along tau.
    double base_re[L], base_im[L], step_re[L] = {}, step_im[L] = {};

    for (size_t j = 0; j < N_tau; ++j) {
      const MathUtils::RealFourVector& u = trajectory[j].momentum;

      // Advance the phasors to tau_j by the small-angle step (branch-free over the block, so it vectorizes),
      // then re-anchor from the exact phase at anchor steps and for any point whose step angle is too large.
      // The step phasor is only needed (and only advanced) when there is more than one frequency.
      bool anchor = j % kAnchorInterval == 0;
      // Harmonic N advances by roughly 2*pi*N/trajectory_NT per step, so high harmonics exceed
      // kMaxSteppedAngle. Then the polynomial is evaluated at angle/2^m and squared m times (same m for the
      // whole block, so the loops stay vectorizable); beyond kMaxHalvings the block falls back to std::polar.
      bool any_large_angle = false;
      if (!anchor) {
        double phi_increment[L];
        double max_increment = 0.0;
        for (size_t l = 0; l < L; ++l) {
          phi_increment[l] = dt0[j] - (nx[l] * drx[j] + ny[l] * dry[j] + nz[l] * drz[j]);
          max_increment = std::max(max_increment, std::abs(phi_increment[l]));
        }
        double max_angle = max_increment * std::max(std::abs(k0), multi_frequency ? std::abs(step) : 0.0);
        size_t halvings = 0;
        double angle_scale = 1.0;
        while (max_angle * angle_scale > kMaxSteppedAngle && halvings < kMaxHalvings) {
          angle_scale *= 0.5;
          ++halvings;
        }
        any_large_angle = max_angle * angle_scale > kMaxSteppedAngle;
        if (!any_large_angle) {
          advance_phasors(base_re, base_im, phi_increment, k0 * angle_scale, halvings);
          if (multi_frequency) advance_phasors(step_re, step_im, phi_increment, step * angle_scale, halvings);
        }
      }
      if (anchor || any_large_angle) {
        for (size_t l = 0; l < L; ++l) {
          double phase = t0[j] - (nx[l] * rx[j] + ny[l] * ry[j] + nz[l] * rz[j]);
          std::complex<double> base = std::polar(1.0, phase * k0);
          base_re[l] = base.real();
          base_im[l] = base.imag();
          if (multi_frequency) {
            std::complex<double> step_phasor = std::polar(1.0, phase * step);
            step_re[l] = step_phasor.real();
            step_im[l] = step_phasor.imag();
          }
        }
      }

      double tau_weight = d_tau;
      if (N_tau > 1 && (j == 0 || j == N_tau - 1)) tau_weight *= 0.5;

      // Frequency-independent factors per screen point (see compute_radiation_long_distance).
      double bulk_term[2][L], endpoint_term[2][L] = {};
      double endpoint_sign = 0.0;
      if (!use_direct_formula) {
        if (j == 0) endpoint_sign += 1.0;
        if (j == N_tau - 1) endpoint_sign -= 1.0;
      }
      bool has_endpoint_contribution = endpoint_sign != 0.0;
      for (size_t l = 0; l < L; ++l) {
        double n0_contract_u = u[0] - (nx[l] * u[1] + ny[l] * u[2] + nz[l] * u[3]);
        double e1_dot_u = e1x[l] * u[1] + e1y[l] * u[2] + e1z[l] * u[3];
        double e2_dot_u = e2x[l] * u[1] + e2y[l] * u[2] + e2z[l] * u[3];
        if (use_direct_formula) {
          const MathUtils::RealFourVector& w = trajectory[j].acceleration;
          double n0_contract_w = w[0] - (nx[l] * w[1] + ny[l] * w[2] + nz[l] * w[3]);
          double e1_dot_w = e1x[l] * w[1] + e1y[l] * w[2] + e1z[l] * w[3];
          double e2_dot_w = e2x[l] * w[1] + e2y[l] * w[2] + e2z[l] * w[3];
          double scale = -tau_weight / (n0_contract_u * n0_contract_u);
          bulk_term[0][l] = scale * (n0_contract_u * e1_dot_w - n0_contract_w * e1_dot_u);
          bulk_term[1][l] = scale * (n0_contract_u * e2_dot_w - n0_contract_w * e2_dot_u);
        } else {
          bulk_term[0][l] = tau_weight * e1_dot_u;
          bulk_term[1][l] = tau_weight * e2_dot_u;
        }
        if (has_endpoint_contribution) {
          endpoint_term[0][l] = endpoint_sign * e1_dot_u / n0_contract_u;
          endpoint_term[1][l] = endpoint_sign * e2_dot_u / n0_contract_u;
        }
      }

      double freq_re[L], freq_im[L];
      for (size_t l = 0; l < L; ++l) {
        freq_re[l] = base_re[l];
        freq_im[l] = base_im[l];
      }
      for (size_t i_freq = 0; i_freq < N_freq; ++i_freq) {
        for (size_t a = 0; a < 2; ++a) {
          size_t offset = (i_freq * 2 + a) * L;
          for (size_t l = 0; l < L; ++l) {
            bulk_re[offset + l] += bulk_term[a][l] * freq_re[l];
            bulk_im[offset + l] += bulk_term[a][l] * freq_im[l];
          }
          if (has_endpoint_contribution) {
            for (size_t l = 0; l < L; ++l) {
              end_re[offset + l] += endpoint_term[a][l] * freq_re[l];
              end_im[offset + l] += endpoint_term[a][l] * freq_im[l];
            }
          }
        }
        for (size_t l = 0; l < L; ++l) {
          double re = freq_re[l] * step_re[l] - freq_im[l] * step_im[l];
          freq_im[l] = freq_re[l] * step_im[l] + freq_im[l] * step_re[l];
          freq_re[l] = re;
        }
      }
    }

    // Reconstruction, identical to compute_radiation_long_distance's.
    for (size_t l = 0; l < block_size; ++l) {
      size_t i_d = block_begin + l;
      const std::array<std::array<double, 3>, 2> basis{{{e1x[l], e1y[l], e1z[l]}, {e2x[l], e2y[l], e2z[l]}}};
      std::array<std::array<double, 6>, 2> basis_bivector{};
      for (size_t a = 0; a < 2; ++a) {
        const std::array<double, 3>& e = basis[a];
        basis_bivector[a] = {
            -e[0], -e[1], -e[2], e[0] * ny[l] - nx[l] * e[1], e[0] * nz[l] - nx[l] * e[2], e[1] * nz[l] - ny[l] * e[2]};
      }
      double inv_x0_norm = 1.0 / D[l];
      double constant_phase = D[l] + reference[0];  // |x_0| + r^0(tau_m), kept out of the stepped phase
      for (size_t i_freq = 0; i_freq < N_freq; i_freq++) {
        MathUtils::Complex amplitude = std::polar(inv_x0_norm, frequencies_list[i_freq] * constant_phase);
        MathUtils::Complex bulk_factor =
            use_direct_formula ? amplitude : MathUtils::Complex{0.0, frequencies_list[i_freq]} * amplitude;
        Simulation::PackedFaraday& target = field.field[i_freq][i_d];
        for (size_t a = 0; a < 2; ++a) {
          size_t index_acc = (i_freq * 2 + a) * L + l;
          MathUtils::Complex p_bulk = bulk_factor * MathUtils::Complex{bulk_re[index_acc], bulk_im[index_acc]};
          MathUtils::Complex p_endpoint = amplitude * MathUtils::Complex{end_re[index_acc], end_im[index_acc]};
          for (size_t index = 0; index < 6; ++index) {
            target.long_range[index] += p_bulk * basis_bivector[a][index];
            target.boundary[index] += p_endpoint * basis_bivector[a][index];
          }
        }
      }
    }
  }
}

}  // namespace

RadiationFormula parse_radiation_formula(const std::string& name) {
  if (name == "simplified") return RadiationFormula::Simplified;
  if (name == "direct") return RadiationFormula::Direct;
  if (name == "long_distance_simplified") return RadiationFormula::LongDistanceSimplified;
  if (name == "long_distance_direct") return RadiationFormula::LongDistanceDirect;
  if (name == "long_distance_simplified_approx") return RadiationFormula::LongDistanceSimplifiedApprox;
  if (name == "long_distance_direct_approx") return RadiationFormula::LongDistanceDirectApprox;
  throw std::runtime_error("unknown radiation_formula \"" + name +
                           "\" (expected simplified, direct, long_distance_simplified, long_distance_direct, "
                           "long_distance_simplified_approx or long_distance_direct_approx)");
}

void compute_radiation(Particle::Electron& electron, const Laser::LaserField& laser,
                       const std::vector<double>& frequencies_list, const Detector::Detector_2D& detector,
                       Simulation::PackedRadiationField& field, RadiationFormula formula) {
  electron.compute_trajectory(laser);
  switch (formula) {
    case RadiationFormula::Simplified:
      compute_radiation_exact(electron, frequencies_list, detector, field, false);
      break;
    case RadiationFormula::Direct:
      compute_radiation_exact(electron, frequencies_list, detector, field, true);
      break;
    case RadiationFormula::LongDistanceSimplified:
      compute_radiation_long_distance(electron, frequencies_list, detector, field, false);
      break;
    case RadiationFormula::LongDistanceDirect:
      compute_radiation_long_distance(electron, frequencies_list, detector, field, true);
      break;
    case RadiationFormula::LongDistanceSimplifiedApprox:
      compute_radiation_long_distance_stepped(electron, frequencies_list, detector, field, false);
      break;
    case RadiationFormula::LongDistanceDirectApprox:
      compute_radiation_long_distance_stepped(electron, frequencies_list, detector, field, true);
      break;
  }
}

}  // namespace Core::Radiation
