#pragma once
#include <complex>
#include <cstddef>
#include <vector>

#include "../detector/detector.hpp"
#include "../io_utils/config_map.hpp"
#include "../laser/laser_field.hpp"
#include "../math_utils/math_utils.hpp"
#include "../particle/electron.hpp"

namespace Core::Simulation {

struct simulation_parameters {
  double tau_0_traj;  // the initial value of the proper time for the electron trajecotry simulation
  double d_tau_traj;
  size_t simulation_length;  // the total number of points in the numerical integration of the electron trajectory
  std::vector<double> frequencies;  // the vector of calculated frequencies.

  // The fundamental's frequency (PhysUtils::non_linear_Thomson_formula(k1, p, n2, 1), i.e. omega/c units,
  // same convention as `frequencies`), computed once regardless of dense_frequency_spectrum -- used by
  // Radiation::plot_radiation_field to normalize the exported "omega" column into units of the fundamental
  // (so a value of 3.0 means "third harmonic"), rather than raw atomic-unit omega.
  double fundamental_frequency;
};

simulation_parameters init_simulation_parameters(const ConfigMap& config, const Laser::LaserField& laser);

// One screen point's contribution at one frequency: the sum, over every electron in the beam, of that
// electron's far-field long-range (velocity) / short-range (acceleration) Faraday tensor F^{\mu\nu} —
// summed as complex amplitudes (not intensities), since capturing the inter-electron interference is
// the whole point of a *coherent* Thomson scattering calculation. Each ComplexFourTensor follows the
// same F^{\mu\nu} packing Laser::LaserField::get_faraday_tensor uses (E = F[0][1..3], B =
// F[1][2]/F[2][3]/F[3][1], antisymmetric elsewhere) — only complex-valued, since these are
// frequency-domain amplitudes rather than real time-domain field values. This is the true,
// fully-reconstructed tensor, built (once, via MathUtils::unpack_bivector) only after every
// electron's and every thread's contribution has been summed — see PackedFaraday for the
// accumulation-time representation. Observables should be computed from this, not PackedFaraday.
struct Faraday {
  Core::MathUtils::ComplexFourTensor long_range;
  Core::MathUtils::ComplexFourTensor short_range;
};

// The coherently-summed radiation field, indexed field[i_omega][i_screen_point].
struct RadiationField {
  std::vector<std::vector<Faraday>> field;
};

// Accumulation-time representation of one screen point's contribution at one frequency: same
// long_range/short_range Faraday split as above, but each held as the packed 6-element
// MathUtils::ComplexBivector rather than the full 4x4 tensor, since compute_radiation's tau loop
// rewrites these at O(N_tau) frequency per screen point and the diagonal/mirrored lower-triangle
// of the corresponding Faraday tensor are never touched until accumulation is finished. Used by
// Radiation::compute_radiation and run_simulation's per-thread accumulators; run_simulation
// unpacks into a Faraday/RadiationField (via MathUtils::unpack_bivector) once reduction across
// threads is complete.
struct PackedFaraday {
  Core::MathUtils::ComplexBivector long_range;
  Core::MathUtils::ComplexBivector short_range;
};

struct PackedRadiationField {
  std::vector<std::vector<PackedFaraday>> field;
};

// Computes the beam's radiated field on the detector screen, over the configured frequency grid.
// Partitions the electron beam across a pool of worker threads; each thread accumulates its assigned
// electrons' contributions into a private RadiationField (avoiding any cross-thread contention during
// the per-electron work), and the thread-local fields are summed once all threads finish.
// num_threads == 0 auto-detects via std::thread::hardware_concurrency(); a value larger than the
// hardware maximum, or larger than the electron beam size, is clamped back down. Pass an explicit
// value (e.g. 1) to force a specific thread count. num_threads is taken by reference and overwritten
// with the actual thread count used, so callers can report the real value rather than the requested one.
RadiationField run_simulation(const ConfigMap& config, const Laser::LaserField& laser,
                              const Detector::Detector_2D& detector, std::vector<Particle::Electron>& electron_beam,
                              const std::vector<double>& frequencies_list, size_t& num_threads);

}  // namespace Core::Simulation