#include "../include/logging/run_log.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <tuple>

#include "../include/io_utils/io_utils.hpp"
#include "../include/math_utils/math_constants.hpp"
#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Logging {

namespace {

// Reads "<key> <value> [unit]" out of the config and returns the value already converted to atomic
// units -- the same pattern detector_factory.cpp/electron_factory.cpp use to build the objects this
// log describes, repeated here (rather than reused from a shared helper) since read_scaled's only
// job is this one conversion, letting each call site stay a single readable line below.
double read_scaled(const ConfigMap& config, const std::string& key) {
  auto [val, unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, key.c_str()));
  return val * IoUtils::convert_unit_to_number(unit, config);
}

// SI-unit conversions (PhysUtils::AtomicUnits::bohr_radius_m/atomic_time_unit_s -- see that file for
// the CODATA basis), used by write_kv's mm column below and directly for the handful of quantities
// (wavelength in nm, period/pulse duration in fs) that need a different SI unit than the generic
// length-in-mm one write_kv provides -- see CLAUDE.md's "Run log" TODO for why the SI column exists.
double length_au_to_mm(double value_au) { return value_au * PhysUtils::AtomicUnits::bohr_radius_m * 1.0e3; }
double length_au_to_nm(double value_au) { return value_au * PhysUtils::AtomicUnits::bohr_radius_m * 1.0e9; }
double time_au_to_fs(double value_au) { return value_au * PhysUtils::AtomicUnits::atomic_time_unit_s * 1.0e15; }

void write_kv(std::ofstream& file, const std::string& label, double value_au, double lambda_au) {
  file << "  " << std::left << std::setw(28) << label << std::right << std::scientific << std::setprecision(6)
       << std::setw(15) << value_au << " a.u.   " << std::setw(15) << (value_au / lambda_au) << " lambda   "
       << std::setw(15) << length_au_to_mm(value_au) << " mm\n";
}

void write_kv_au(std::ofstream& file, const std::string& label, double value_au) {
  file << "  " << std::left << std::setw(28) << label << std::right << std::scientific << std::setprecision(6)
       << std::setw(15) << value_au << " a.u.\n";
}

void write_section(std::ofstream& file, const std::string& title) {
  file << "\n" << title << "\n" << std::string(title.size(), '-') << "\n";
}

}  // namespace

void write_run_log(const ConfigMap& config, const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                   size_t num_electrons, const Simulation::simulation_parameters& sim_par, size_t num_threads,
                   double simulation_elapsed_seconds, double total_cpu_seconds, const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for run log export: " + filepath);
  }

  const double c = PhysUtils::AtomicUnits::c;
  const double mc = PhysUtils::AtomicUnits::m_0 * c;
  const double lambda_au = 2.0 * MathUtils::pi / laser.get_omega() * c;
  const std::string laser_type = IoUtils::get_required(config, "laser_type");
  const bool is_laguerre_gauss = (laser_type == "laguerre_gauss");
  double lg_w0_au = 0.0;
  if (is_laguerre_gauss) {
    lg_w0_au = std::get<2>(IoUtils::get_laser_lg_params(config));
  }

  file << "CoherentThomson run log\n"
       << "========================\n"
       << "Everything here is derived from config.cfg (also copied verbatim into this run directory) plus what "
          "the solver actually computed from it -- meant for identifying at a glance what a given run's "
          "output actually corresponds to, without re-deriving units/derived quantities by hand.\n";

  write_section(file, "Laser");
  file << "  type                        " << laser_type << "\n";
  write_kv_au(file, "omega", laser.get_omega());
  write_kv_au(file, "lambda", lambda_au);
  file << "  lambda (nm)                 " << length_au_to_nm(lambda_au) << " nm\n";
  double period_T_au = IoUtils::get_laser_period(config);
  write_kv_au(file, "period T", period_T_au);
  file << "  period T (fs)               " << time_au_to_fs(period_T_au) << " fs\n";
  file << "  a0                          " << laser.get_a0() << "\n";
  file << "  zeta_1                      (" << laser.get_zeta_1().real() << ", " << laser.get_zeta_1().imag() << ")\n";
  file << "  zeta_2                      (" << laser.get_zeta_2().real() << ", " << laser.get_zeta_2().imag() << ")\n";
  MathUtils::RealFourVector laser_dir = IoUtils::get_laser_direction(config);
  file << "  direction (nx,ny,nz)        (" << laser_dir[1] << ", " << laser_dir[2] << ", " << laser_dir[3]
       << ")  [raw, pre-normalization]\n";
  write_kv_au(file, "flat_duration (phase)", laser.get_flat_duration());
  file << "  flat_duration (cycles)      " << laser.get_flat_duration() / (2.0 * MathUtils::pi) << "\n";
  write_kv_au(file, "wing_sigma (phase)", IoUtils::get_laser_wing_sigma(config));
  file << "  wing_sigma_cutoff           " << laser.get_wing_sigma_cutoff() << "\n";
  // Total pulse duration (flat-top + both wings, sim_par.total_cycles -- the same quantity used to
  // size the electron trajectory grid in Simulation::init_simulation_parameters and, further
  // normalized, the "CPU time /electron/screen pt/cycle" figure below), unlike flat_duration above
  // which only covers the flat-top part.
  double pulse_duration_au = sim_par.total_cycles * period_T_au;
  file << "  pulse duration (cycles)     " << sim_par.total_cycles << "\n";
  write_kv_au(file, "pulse duration", pulse_duration_au);
  file << "  pulse duration (fs)         " << time_au_to_fs(pulse_duration_au) << " fs\n";
  if (is_laguerre_gauss) {
    auto [lg_p, lg_l, lg_w0] = IoUtils::get_laser_lg_params(config);
    file << "  LG p, l                     " << lg_p << ", " << lg_l << "\n";
    write_kv(file, "LG w0", lg_w0, lambda_au);
  }

  write_section(file, "Electron beam");
  file << "  particle count              " << num_electrons << "\n";
  IoUtils::CylinderBeamParams beam = IoUtils::parse_cylinder_beam_params(config);
  write_kv(file, "cylinder radius", beam.radius, lambda_au);
  write_kv(file, "cylinder height", beam.height, lambda_au);
  write_kv(file, "center_x", beam.center_x, lambda_au);
  write_kv(file, "center_y", beam.center_y, lambda_au);
  write_kv(file, "center_z", beam.center_z, lambda_au);
  file << "  average_p (px,py,pz)        (" << beam.average_px / mc << ", " << beam.average_py / mc << ", "
       << beam.average_pz / mc << ")  mc\n";
  file << "  sigma_p (px,py,pz)          (" << beam.sigma_px / mc << ", " << beam.sigma_py / mc << ", "
       << beam.sigma_pz / mc << ")  mc\n";

  write_section(file, "Detector");
  std::string detector_type = detector.get_type_name();
  file << "  type                        " << detector_type << "\n";
  file << "  grid                        " << detector.get_cols() << " x " << detector.get_rows() << "  ("
       << detector.get_total_points() << " points)\n";
  auto [dir_theta, dir_phi] = IoUtils::get_detector_direction_angles(config);
  file << "  direction (theta,phi)       (" << dir_theta / MathUtils::pi << ", " << dir_phi / MathUtils::pi
       << ")  pi\n";
  if (detector_type == "RectangularDetector") {
    write_kv(file, "distance", read_scaled(config, "rectangular_detector_distance"), lambda_au);
    write_kv(file, "x_min", read_scaled(config, "rectangular_detector_x_min"), lambda_au);
    write_kv(file, "x_max", read_scaled(config, "rectangular_detector_x_max"), lambda_au);
    write_kv(file, "y_min", read_scaled(config, "rectangular_detector_y_min"), lambda_au);
    write_kv(file, "y_max", read_scaled(config, "rectangular_detector_y_max"), lambda_au);
    if (is_laguerre_gauss) {
      file << "  x_max / w0                  " << read_scaled(config, "rectangular_detector_x_max") / lg_w0_au << "\n";
      file << "  y_max / w0                  " << read_scaled(config, "rectangular_detector_y_max") / lg_w0_au << "\n";
    }
  } else if (detector_type == "SphericalDetector") {
    write_kv(file, "radius", read_scaled(config, "spherical_detector_radius"), lambda_au);
    file << "  theta_min, theta_max        (" << read_scaled(config, "spherical_detector_theta_min") / MathUtils::pi
         << ", " << read_scaled(config, "spherical_detector_theta_max") / MathUtils::pi << ")  pi\n";
    file << "  phi_min, phi_max            (" << read_scaled(config, "spherical_detector_phi_min") / MathUtils::pi
         << ", " << read_scaled(config, "spherical_detector_phi_max") / MathUtils::pi << ")  pi\n";
  } else if (detector_type == "CircularDetector") {
    write_kv(file, "distance", read_scaled(config, "circular_detector_distance"), lambda_au);
    write_kv(file, "R_min", read_scaled(config, "circular_detector_R_min"), lambda_au);
    write_kv(file, "R_max", read_scaled(config, "circular_detector_R_max"), lambda_au);
    if (is_laguerre_gauss) {
      file << "  R_max / w0                  " << read_scaled(config, "circular_detector_R_max") / lg_w0_au << "\n";
    }
  }

  write_section(file, "Frequency spectrum");
  const std::string& radiation_formula = IoUtils::get_required(config, "radiation_formula");
  file << "  radiation_formula           " << radiation_formula << "\n";
  bool dense_spectrum = IoUtils::get_required(config, "dense_frequency_spectrum") == "true";
  file << "  dense_frequency_spectrum    " << (dense_spectrum ? "true" : "false") << "\n";
  if (dense_spectrum) {
    IoUtils::OmegaRange omega_range = IoUtils::get_omega_range(config);
    auto unit_name = [](IoUtils::OmegaRangeUnit u) {
      return u == IoUtils::OmegaRangeUnit::FirstHarmonicFrequency ? "first_harmonic_frequency" : "omega_laser";
    };
    file << "  N_omega                     " << IoUtils::get_number_of_frequencies(config) << "\n";
    file << "  omega_min (as configured)   " << omega_range.omega_min << "  [" << unit_name(omega_range.min_unit)
         << "]\n";
    file << "  omega_max (as configured)   " << omega_range.omega_max << "  [" << unit_name(omega_range.max_unit)
         << "]\n";
  } else {
    file << "  N_harmonics                 " << IoUtils::get_number_of_harmonics(config) << "\n";
    file << "  N_harmonics_min             " << IoUtils::get_number_of_harmonics_min(config) << "\n";
  }
  write_kv_au(file, "fundamental_frequency", sim_par.fundamental_frequency * c);
  file << "  fundamental / laser omega   " << sim_par.fundamental_frequency * c / laser.get_omega()
       << "  [Doppler/nonlinear shift ratio]\n";
  file << "  dressed momentum q          (" << sim_par.q[0] << ", " << sim_par.q[1] << ", " << sim_par.q[2] << ", "
       << sim_par.q[3] << ")  a.u.\n";

  file << "  frequencies_list (" << sim_par.frequencies.size() << " points, omega/c convention, "
       << "normalized by fundamental_frequency to match radiation_field.dat's own omega column):\n";
  if (sim_par.frequencies.size() <= 20) {
    for (size_t i = 0; i < sim_par.frequencies.size(); ++i) {
      file << "    [" << i << "] " << sim_par.frequencies[i] << "  ("
           << sim_par.frequencies[i] / sim_par.fundamental_frequency << " omega_1)\n";
    }
  } else {
    double first = sim_par.frequencies.front();
    double last = sim_par.frequencies.back();
    double step =
        (sim_par.frequencies.size() > 1) ? (last - first) / static_cast<double>(sim_par.frequencies.size() - 1) : 0.0;
    file << "    first = " << first << "  (" << first / sim_par.fundamental_frequency << " omega_1)\n";
    file << "    last  = " << last << "  (" << last / sim_par.fundamental_frequency << " omega_1)\n";
    file << "    step  = " << step << "  (evenly spaced)\n";
  }
  file << "  print_field_in_canonical_frame  " << IoUtils::get_required(config, "print_field_in_canonical_frame")
       << "\n";

  write_section(file, "Run");
  file << "  num_threads (resolved)      " << num_threads << "\n";
  file << "  random_seed                 " << IoUtils::get_required(config, "random_seed") << "\n";
  file << "  debug                       " << IoUtils::get_required(config, "debug") << "\n";
  file << "  screen points               " << detector.get_total_points() << "\n";
  file << "  simulation wall time (s)    " << simulation_elapsed_seconds << "\n";

  // Per-unit-of-work CPU cost, further normalized by the pulse's own total duration in laser cycles
  // (sim_par.total_cycles) on top of the electron/screen-point normalization run_simulation already
  // prints to stdout -- lets a run's cost be extrapolated to a longer/shorter pulse (more/fewer
  // cycles) directly, rather than just electron count/screen resolution. Built from total_cpu_seconds
  // (summed across every thread's own elapsed time), not simulation_elapsed_seconds -- see this
  // function's doc comment for why. Expressed in microseconds (not seconds) since dividing all the
  // way down to a single electron/screen-point/cycle makes the seconds figure inconveniently small.
  size_t N_screen = detector.get_total_points();
  double time_per_electron_screen_pt_cycle_us = 0.0;
  if (num_electrons > 0 && N_screen > 0 && sim_par.total_cycles > 0.0) {
    double time_per_electron_screen_pt_s =
        total_cpu_seconds / static_cast<double>(num_electrons) / static_cast<double>(N_screen);
    time_per_electron_screen_pt_cycle_us = time_per_electron_screen_pt_s / sim_par.total_cycles * 1.0e6;
  }
  file << "  CPU time /electron/screen pt/cycle   " << time_per_electron_screen_pt_cycle_us << " usec\n";
  file << "  number of frequencies                " << sim_par.frequencies.size() << "\n";
}

}  // namespace Core::Logging
