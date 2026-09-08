
#include "../include/simulation/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

#include "../include/io_utils/io_utils.hpp"
#include "../include/math_utils/math_constants.hpp"
#include "../include/phys_utils/phys_utils.hpp"
#include "../include/radiation/radiation.hpp"

namespace Core::Simulation {

// ComplexBivector aliases std::array, whose only associated namespace (for ADL) is std -- pull
// MathUtils::operator+= into scope explicitly so `bivector += bivector` below resolves.
using Core::MathUtils::operator+=;

simulation_parameters init_simulation_parameters(const ConfigMap& config, const Laser::LaserField& laser) {
  double tau_0_traj = 0.0;  // hardcoded, to be modified and included in CLAUDE.md
  double T = IoUtils::get_laser_period(config);
  size_t NT_traj = IoUtils::get_trajectory_NT(config);
  double d_tau_traj = T / NT_traj;

  // Total number of proper-time points for the electron trajectory integration: NT_traj points per
  // cycle, times the total number of cycles spanned by the laser pulse. total_phase is read straight
  // off the already-built laser (phi_max - phi_min) rather than re-derived from flat_duration/
  // wing_sigma/wing_sigma_cutoff, so this stays correct even if LaserField's own phase-range formula
  // changes.
  double total_phase = laser.get_phi_max() - laser.get_phi_min();
  double total_cycles = total_phase / (2 * MathUtils::pi);
  size_t simulation_length = static_cast<size_t>(std::round(NT_traj * total_cycles));

  // here we scale the time_step d_tau to accomodate electrons initially in motion.
  // the scale factor is (from the Doppler transofrmation of frequencies)
  // use the formula  d_tau /= (p^0-p^3)/mc where p^0, p^3 are the components of (average) electron intial four-momentum
  IoUtils::CylinderBeamParams beam_params = IoUtils::parse_cylinder_beam_params(config);
  double px = beam_params.average_px;
  double py = beam_params.average_py;
  double pz = beam_params.average_pz;
  double mc = PhysUtils::AtomicUnits::m_0 * PhysUtils::AtomicUnits::c;
  double p0 = std::sqrt(px * px + py * py + pz * pz + mc * mc);
  d_tau_traj /= (p0 - pz) / mc;

  // N_frequencies is read from a different config key depending on the mode: N_harmonics (default,
  // dense_frequency_spectrum=false) vs. N_omega (dense_frequency_spectrum=true) -- see
  // IoUtils::get_number_of_harmonics/get_number_of_frequencies. Kept as two separate keys so
  // switching modes doesn't silently reinterpret whatever count was already configured for the
  // other mode.
  bool dense_spectrum = IoUtils::get_required(config, "dense_frequency_spectrum") == "true";
  size_t N_frequencies =
      dense_spectrum ? IoUtils::get_number_of_frequencies(config) : IoUtils::get_number_of_harmonics(config);
  std::vector<double> frequencies_list(N_frequencies);

  MathUtils::RealFourVector p(p0, px, py, pz);

  // k1, p, and n2 are all evaluated in the canonical frame (laser along Oz), rather than the true lab
  // frame: since the beam and detector are rotated into the lab frame by the same shared rotation as the
  // laser, and non_linear_Thomson_formula only ever combines these through Minkowski contractions (which
  // are invariant under a common rotation of all their arguments), the result is identical either way --
  // and the canonical frame needs no rotation at all. k1 is therefore fixed along canonical Oz (theta=0,
  // phi=0) rather than laser.get_unity_n() (the laser's actual, rotated direction), and n2 is the
  // detector's own canonical-frame direction (detector_direction_theta/phi) rather than the electron's
  // own direction of motion. The actual spectrum depends on the observation direction, but in the code we
  // calculate it at the same set of values for the entire screen. Hard coded to be N_harmonics
  // consecutive harmonics starting at N_harmonics_min by default; the omega_min/omega_max limits in
  // the input file are only honored when dense_frequency_spectrum is true (see CLAUDE.md).
  MathUtils::RealFourVector k1 =
      MathUtils::create_unit_light_like_vector<double>(0.0, 0.0) * laser.get_omega() / PhysUtils::AtomicUnits::c;
  auto [detector_dir_theta, detector_dir_phi] = IoUtils::get_detector_direction_angles(config);
  MathUtils::RealFourVector n2 = MathUtils::create_unit_light_like_vector<double>(detector_dir_theta, detector_dir_phi);

  // compute the dressed (ponderomotive drift) electron momentum q = p + (mc)^2 xi^2 / (2 k1.p) * k1,
  // giving the mass-shell shift m_eff^2 = m^2(1+xi^2/2). The (mc)^2 factor is essential: xi is
  // dimensionless, so without it the correction term has the wrong units and (since mc = 137.036 in
  // atomic units) ends up ~(mc)^2 too small to have any visible effect.
  double xi = laser.get_a0();
  MathUtils::FourVector q = p + mc * mc * xi * xi / (2.0 * MathUtils::contract(p, k1)) * k1;

  // The fundamental's frequency, computed the same way regardless of dense_spectrum -- used only to
  // let Radiation::plot_radiation_field normalize its exported "omega" column into units of the
  // fundamental (see simulation_parameters::fundamental_frequency).
  // to include non linear effects the non_linear_Thomson_formula is called with q
  double fundamental_frequency = PhysUtils::non_linear_Thomson_formula(k1, q, n2, 1);

  // Default: N_harmonics consecutive harmonics starting at N_harmonics_min (1, the common case,
  // reproduces the previous "first N_harmonics harmonics" behavior), for imaging over the whole
  // detector screen. When dense_frequency_spectrum is true, frequencies_list is instead a fine
  // linear scan from omega_min to omega_max (N_omega points) -- meant for a detector collapsed to
  // a single point (N_total_points == 1, see main.cpp's warning below), to resolve a Thomson
  // line's width rather than just locate the harmonic peaks.
  if (dense_spectrum) {
    // the scaling factor has the role of adjusting the frequency interval position
    // in the case when the emitted frequency is not multiple of omega_laser (this includes also non linear effects)
    double frequency_scaling_factor = fundamental_frequency * PhysUtils::AtomicUnits::c / laser.get_omega();
    auto [omega_min, omega_max] = IoUtils::get_omega_range(config);
    omega_min *= frequency_scaling_factor;
    omega_max *= frequency_scaling_factor;
    std::cout << "frequency_scaling_factor " << frequency_scaling_factor << std::endl;
    for (size_t i = 0; i < N_frequencies; i++) {
      double omega = (N_frequencies > 1) ? omega_min + static_cast<double>(i) * (omega_max - omega_min) /
                                                           static_cast<double>(N_frequencies - 1)
                                         : omega_min;
      frequencies_list[i] = omega / PhysUtils::AtomicUnits::c;  // match non_linear_Thomson_formula's omega/c convention
    }
  } else {
    size_t N_harmonics_min = IoUtils::get_number_of_harmonics_min(config);
    for (size_t i = 0; i < N_frequencies; i++) {
      frequencies_list[i] = PhysUtils::non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i);
    }
  }

  return simulation_parameters{tau_0_traj, d_tau_traj, simulation_length, frequencies_list, fundamental_frequency};
};

RadiationField run_simulation(const ConfigMap& config, const Laser::LaserField& laser,
                              const Detector::Detector_2D& detector, std::vector<Particle::Electron>& electron_beam,
                              const std::vector<double>& frequencies_list, size_t& num_threads) {
  size_t N_omega = frequencies_list.size();
  size_t N_screen = detector.get_total_points();
  size_t num_electrons = electron_beam.size();

  // Selects which of the two closed forms in theory/FT_Faraday_tensor-direct_and_simplified_forms.md
  // Radiation::compute_radiation evaluates -- see that function's doc comment. Validated here
  // (rather than left as a silent `== "direct"` comparison) so a typo'd value doesn't silently fall
  // back to "simplified" the way an unrecognized config unit silently falls back to 1.0 elsewhere.
  const std::string& radiation_formula = IoUtils::get_required(config, "radiation_formula");
  if (radiation_formula != "simplified" && radiation_formula != "direct") {
    throw std::runtime_error("Simulation::run_simulation: unknown radiation_formula \"" + radiation_formula +
                             "\" (expected \"simplified\" or \"direct\")");
  }
  bool use_direct_formula = radiation_formula == "direct";

  size_t max_threads = std::thread::hardware_concurrency();
  if (max_threads == 0) max_threads = 1;
  if (num_threads == 0 || num_threads > max_threads) num_threads = max_threads;
  num_threads = std::min(num_threads, num_electrons);

  // One private accumulator per thread, each grid sized [N_omega][N_screen] and zero-initialized,
  // holding the packed (6-complex-element) accumulation-time Faraday representation (see
  // PackedFaraday) so each thread can accumulate its assigned electrons' contributions without
  // touching any other thread's memory, and without paying for the full 4x4 tensor's
  // storage/bandwidth until the true Faraday tensor is built below.
  std::vector<PackedRadiationField> thread_fields(
      num_threads,
      PackedRadiationField{std::vector<std::vector<PackedFaraday>>(N_omega, std::vector<PackedFaraday>(N_screen))});

  // Split the beam into num_threads contiguous chunks of electrons, one chunk per thread.
  size_t chunk_size = (num_electrons + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t begin = thread_idx * chunk_size;
    size_t end = std::min(begin + chunk_size, num_electrons);
    if (begin >= end) continue;

    threads.emplace_back([&, begin, end, thread_idx]() {
      PackedRadiationField& local_field = thread_fields[thread_idx];
      for (size_t p = begin; p < end; ++p) {
        Particle::Electron& electron = electron_beam[p];
        Radiation::compute_radiation(electron, laser, frequencies_list, detector, local_field, use_direct_formula);
        // Progress indicator: only thread 0 prints, both to avoid interleaved output from multiple
        // threads writing to std::cout concurrently and because thread 0's chunk is representative
        // enough of overall progress for a rough sense of how a long run is advancing.
        if (thread_idx == 0) {
          std::cout << "\rThread 0: electron " << (p - begin + 1) << "/" << (end - begin) << std::flush;
        }
      }
      if (thread_idx == 0) {
        std::cout << "\n";
      }
    });
  }
  for (auto& thread : threads) {
    thread.join();
  }

  // Reduce: sum the thread-local packed fields (amplitudes add linearly, so order doesn't matter)
  // into a combined packed result.
  PackedRadiationField packed_result{
      std::vector<std::vector<PackedFaraday>>(N_omega, std::vector<PackedFaraday>(N_screen))};
  for (const PackedRadiationField& local_field : thread_fields) {
    for (size_t i_omega = 0; i_omega < N_omega; ++i_omega) {
      for (size_t i_screen = 0; i_screen < N_screen; ++i_screen) {
        packed_result.field[i_omega][i_screen].long_range += local_field.field[i_omega][i_screen].long_range;
        packed_result.field[i_omega][i_screen].short_range += local_field.field[i_omega][i_screen].short_range;
      }
    }
  }

  // Build the true Faraday tensor (Simulation::Faraday, full 4x4 F^{mu nu}) from the packed
  // accumulation above -- this is the only place the antisymmetric lower triangle and zero
  // diagonal are ever filled in, once per screen point/frequency (via MathUtils::unpack_bivector)
  // instead of on every trajectory point of every electron. Everything downstream (plotting,
  // observable calculations) consumes this, not the packed representation.
  //
  // compute_radiation builds E/B (and hence F) straight in the lab frame, using n_R0/u vectors that
  // already carry the laser's actual (rotated) orientation -- but plots are usually meant to be read in
  // the laser's own canonical frame (laser along Oz), the same convention init_simulation_parameters
  // uses for k1/p/n2 above. Rotating back with the inverse of the laser's rotation_matrix
  // (MathUtils::inverse_rotation_tensor + rotate_tensor) once here, rather than per trajectory point, is
  // gated behind the "print_field_in_canonical_frame" config key so the lab-frame field (as actually
  // computed, matching laser_nx/ny/nz) can still be inspected when that's what's wanted.
  bool print_field_in_canonical_frame = IoUtils::get_required(config, "print_field_in_canonical_frame") == "true";
  MathUtils::RealFourTensor inverse_rotation = MathUtils::inverse_rotation_tensor(laser.get_rotation_matrix());
  RadiationField result{std::vector<std::vector<Faraday>>(N_omega, std::vector<Faraday>(N_screen))};
  for (size_t i_omega = 0; i_omega < N_omega; ++i_omega) {
    for (size_t i_screen = 0; i_screen < N_screen; ++i_screen) {
      const PackedFaraday& packed = packed_result.field[i_omega][i_screen];
      MathUtils::ComplexFourTensor long_range = MathUtils::unpack_bivector(packed.long_range);
      MathUtils::ComplexFourTensor short_range = MathUtils::unpack_bivector(packed.short_range);
      if (print_field_in_canonical_frame) {
        long_range = MathUtils::rotate_tensor(inverse_rotation, long_range);
        short_range = MathUtils::rotate_tensor(inverse_rotation, short_range);
      }
      result.field[i_omega][i_screen].long_range = long_range;
      result.field[i_omega][i_screen].short_range = short_range;
    }
  }

  // multiply with a common factor; it reduces to 1/(2pi c^2) in atomic units. The c^2 (not c) is
  // required by the Jacobian of the t->tau change of variable (dt/d(tau) = (u.n_R0)/c, on top of
  // Jackson's own 1/c) -- see theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Jacobian
  // of the change of variable" note; both the simplified and direct formulas now share this same
  // overall constant.

  double general_factor =
      1 / (2 * MathUtils::pi) * PhysUtils::AtomicUnits::e_0 /
      (4 * MathUtils::pi * PhysUtils::AtomicUnits::epsilon_0 * PhysUtils::AtomicUnits::c * PhysUtils::AtomicUnits::c);

  for (size_t i_omega = 0; i_omega < N_omega; ++i_omega) {
    for (size_t i_screen = 0; i_screen < N_screen; ++i_screen) {
      result.field[i_omega][i_screen].long_range *= general_factor;
      result.field[i_omega][i_screen].short_range *= general_factor;
    }
  }

  return result;
}

}  // namespace Core::Simulation