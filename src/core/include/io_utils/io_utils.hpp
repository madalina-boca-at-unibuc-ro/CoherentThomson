#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "../math_utils/math_constants.hpp"
#include "../math_utils/math_utils.hpp"
#include "../phys_utils/phys_utils.hpp"
#include "config_map.hpp"

namespace Core::IoUtils {

// Looks up `key` in `config`, throwing a std::out_of_range that names the missing key rather than
// std::map::at's unhelpful "map::at". Use this instead of config.at(key) everywhere a config key is
// read, so a missing/misspelled key in the .cfg file is diagnosable from the error message alone.
inline const std::string& get_required(const ConfigMap& config, const char* key) {
  auto it = config.find(key);
  if (it == config.end()) {
    throw std::out_of_range(std::string("Missing required config key: \"") + key + "\"");
  }
  return it->second;
}

// Helper helper to trim trailing spaces and carriage returns (\r from Windows
// files)
inline std::string trim_trailing(const std::string& str) {
  auto end = std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) { return !std::isspace(ch); });
  return (end == str.rend()) ? "" : std::string(str.begin(), end.base());
}

// Splits a string into its numerical magnitude and its trailing unit identifier
// text
inline std::pair<double, std::string> split_value_and_unit(const std::string& raw_value) {
  std::stringstream ss(raw_value);
  double numeric_value;
  std::string unit_text = "";

  // 1. Extract the number component
  if (!(ss >> numeric_value)) {
    // If it's not a valid number at all, fallback safely to zero or handle it
    return {0.0, "INVALID"};
  }

  // 2. Extract the trailing text string if it exists
  ss >> unit_text;

  // Returns a pair: {magnitude, "unit"}
  return {numeric_value, unit_text};
}

inline std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

// Case-insensitive check of whether unit_name appears in allowed_units.
inline bool is_unit_allowed(const std::string& unit_name, const std::vector<std::string>& allowed_units) {
  std::string lower_unit_name = to_lower(unit_name);
  return std::any_of(allowed_units.begin(), allowed_units.end(),
                     [&](const std::string& allowed) { return to_lower(allowed) == lower_unit_name; });
}

inline std::string convert_bool_to_string(const bool& value) {
  std::string s = value ? "true" : "false";
  return s;
}

inline double convert_unit_to_number(std::string unit_name, const ConfigMap& config,
                                     const std::vector<std::string>& allowed_units = {}) {
  if (!allowed_units.empty() && !is_unit_allowed(unit_name, allowed_units)) {
    throw std::invalid_argument("Unit \"" + unit_name + "\" is not among the allowed units for this value.");
  }

  if (to_lower(unit_name) == "pi") {
    return MathUtils::pi;
  } else if (to_lower(unit_name) == "lambda") {
    double omega = std::stod(get_required(config, "laser_frequency"));
    return 2 * MathUtils::pi / omega * PhysUtils::AtomicUnits::c;
  } else if (to_lower(unit_name) == "mc") {
    return PhysUtils::AtomicUnits::m_0 * PhysUtils::AtomicUnits::c;
  } else if (to_lower(unit_name) == "cycles_adim") {
    return 2 * MathUtils::pi;
  } else if (to_lower(unit_name) == "omega_laser") {
    // Read directly rather than via get_laser_frequency (defined later in this file) to avoid a
    // forward reference -- same pattern the "lambda" branch above already uses.
    return std::stod(get_required(config, "laser_frequency"));
  } else if (to_lower(unit_name) == "a.u.") {
    return 1.0;
  } else {
    std::cerr << "Warning: Unknown unit " << unit_name << ". Assuming 1.0." << std::endl;
    return 1.0;
  }
}

inline double get_laser_frequency(const ConfigMap& config) {
  return std::stod(get_required(config, "laser_frequency"));
}

inline double get_laser_a0(const ConfigMap& config) { return std::stod(get_required(config, "laser_a0")); }

// Half-width of the flat-top plateau, in radians.
inline double get_laser_flat_duration(const ConfigMap& config) {
  auto [value, unit] = split_value_and_unit(get_required(config, "laser_flat_duration"));
  return value * convert_unit_to_number(unit, config);
}

// Standard deviation (phase width) of the Gaussian wings, in radians.
inline double get_laser_wing_sigma(const ConfigMap& config) {
  auto [value, unit] = split_value_and_unit(get_required(config, "laser_wing_sigma"));
  return value * convert_unit_to_number(unit, config);
}

// Global phase shift of the pulse, in radians.
inline double get_laser_delay(const ConfigMap& config) {
  auto [value, unit] = split_value_and_unit(get_required(config, "laser_delay"));
  return value * convert_unit_to_number(unit, config);
}

// Number of Gaussian-wing standard deviations included in the simulated phase range.
inline double get_laser_wing_sigma_cutoff(const ConfigMap& config) {
  return std::stod(get_required(config, "laser_wing_sigma_cutoff"));
}

// Raw (unnormalized) propagation direction, as a 4-vector with a zero time component.
inline MathUtils::RealFourVector get_laser_direction(const ConfigMap& config) {
  double nx = std::stod(get_required(config, "laser_nx"));
  double ny = std::stod(get_required(config, "laser_ny"));
  double nz = std::stod(get_required(config, "laser_nz"));
  return MathUtils::RealFourVector(0.0, nx, ny, nz);
}

// Reads a complex-valued config entry stored as two keys, "<base_key>_re" and "<base_key>_im".
inline MathUtils::Complex get_complex_config_value(const ConfigMap& config, const std::string& base_key) {
  double re = std::stod(get_required(config, (base_key + "_re").c_str()));
  double im = std::stod(get_required(config, (base_key + "_im").c_str()));
  return MathUtils::Complex(re, im);
}

// Complex polarization coefficients for the two orthogonal polarization directions: {zeta_1, zeta_2}.
// Circular polarization is zeta_1 = (1,0), zeta_2 = (0,1); linear along epsilon_1 is zeta_1 = (1,0),
// zeta_2 = (0,0).
inline std::pair<MathUtils::Complex, MathUtils::Complex> get_laser_zeta(const ConfigMap& config) {
  return {get_complex_config_value(config, "laser_zeta_1"), get_complex_config_value(config, "laser_zeta_2")};
}

// Laguerre-Gauss mode parameters: {radial index p, azimuthal/OAM index l, beam waist w0}. Only read
// when laser_type == "laguerre_gauss".
inline std::tuple<int, int, double> get_laser_lg_params(const ConfigMap& config) {
  int p = std::stoi(get_required(config, "laser_lg_p"));
  int l = std::stoi(get_required(config, "laser_lg_l"));
  auto [w0_val, w0_unit] = split_value_and_unit(get_required(config, "laser_lg_w0"));
  double w0 = w0_val * convert_unit_to_number(w0_unit, config);
  return {p, l, w0};
}

// Detector's own normal direction, as {theta, phi} in radians, given in the canonical frame where the
// laser propagates along Oz (independent of laser_nx/ny/nz) — shared by Detector::create_detector and
// Simulation::init_simulation_parameters so both read the same direction.
inline std::pair<double, double> get_detector_direction_angles(const ConfigMap& config) {
  auto [theta_val, theta_unit] = split_value_and_unit(get_required(config, "detector_direction_theta"));
  double theta = theta_val * convert_unit_to_number(theta_unit, config);
  auto [phi_val, phi_unit] = split_value_and_unit(get_required(config, "detector_direction_phi"));
  double phi = phi_val * convert_unit_to_number(phi_unit, config);
  return {theta, phi};
}

// omega_min/omega_max, in raw omega (not omega/c) -- only consulted when dense_frequency_spectrum
// is true (Simulation::init_simulation_parameters), to build a fine linear frequency scan instead
// of the default first-N_omega-harmonics list.
inline std::pair<double, double> get_omega_range(const ConfigMap& config) {
  auto [min_val, min_unit] = split_value_and_unit(get_required(config, "omega_min"));
  double omega_min = min_val * convert_unit_to_number(min_unit, config);
  auto [max_val, max_unit] = split_value_and_unit(get_required(config, "omega_max"));
  double omega_max = max_val * convert_unit_to_number(max_unit, config);
  return {omega_min, omega_max};
}

inline size_t get_laser_NT(const ConfigMap& config) { return std::stoul(get_required(config, "laser_NT")); }

inline size_t get_trajectory_NT(const ConfigMap& config) { return std::stoul(get_required(config, "trajectory_NT")); }

// Number of points in the dense_frequency_spectrum=true fine frequency scan (see get_omega_range) --
// distinct from N_harmonics (get_number_of_harmonics), which counts harmonics in the default mode.
inline size_t get_N_omega(const ConfigMap& config) { return std::stoul(get_required(config, "N_omega")); }

// The per-electron parameters a cylinder beam is generated from: cylinder geometry and Gaussian
// momentum distribution. Deliberately excludes beam_particle_count/random_seed (which control the
// generation loop rather than any individual electron) and store_trajectory (a per-run toggle read
// separately). Shared by Particle::generate_cylinder_beam and Simulation::init_simulation_parameters
// (the latter only needs average_px/py/pz, for the initial-motion Doppler scaling of d_tau).
struct CylinderBeamParams {
  double radius;
  double height;
  double center_x, center_y, center_z;
  double average_px, average_py, average_pz;
  double sigma_px, sigma_py, sigma_pz;
};

// Picks the cylinder-beam keys (with unit conversion) out of the config into a CylinderBeamParams.
inline CylinderBeamParams parse_cylinder_beam_params(const ConfigMap& config) {
  CylinderBeamParams params;

  auto [radius_val, radius_unit] = split_value_and_unit(get_required(config, "beam_cylinder_radius"));
  params.radius = radius_val * convert_unit_to_number(radius_unit, config);
  auto [height_val, height_unit] = split_value_and_unit(get_required(config, "beam_cylinder_height"));
  params.height = height_val * convert_unit_to_number(height_unit, config);

  auto [center_x_val, center_x_unit] = split_value_and_unit(get_required(config, "beam_center_x"));
  params.center_x = center_x_val * convert_unit_to_number(center_x_unit, config);
  auto [center_y_val, center_y_unit] = split_value_and_unit(get_required(config, "beam_center_y"));
  params.center_y = center_y_val * convert_unit_to_number(center_y_unit, config);
  auto [center_z_val, center_z_unit] = split_value_and_unit(get_required(config, "beam_center_z"));
  params.center_z = center_z_val * convert_unit_to_number(center_z_unit, config);

  auto [average_px, average_px_units] = split_value_and_unit(get_required(config, "average_px"));
  params.average_px = average_px * convert_unit_to_number(average_px_units, config);
  auto [average_py, average_py_units] = split_value_and_unit(get_required(config, "average_py"));
  params.average_py = average_py * convert_unit_to_number(average_py_units, config);
  auto [average_pz, average_pz_units] = split_value_and_unit(get_required(config, "average_pz"));
  params.average_pz = average_pz * convert_unit_to_number(average_pz_units, config);

  auto [sigma_px, sigma_px_units] = split_value_and_unit(get_required(config, "sigma_px"));
  params.sigma_px = sigma_px * convert_unit_to_number(sigma_px_units, config);
  auto [sigma_py, sigma_py_units] = split_value_and_unit(get_required(config, "sigma_py"));
  params.sigma_py = sigma_py * convert_unit_to_number(sigma_py_units, config);
  auto [sigma_pz, sigma_pz_units] = split_value_and_unit(get_required(config, "sigma_pz"));
  params.sigma_pz = sigma_pz * convert_unit_to_number(sigma_pz_units, config);

  return params;
}

inline double get_lambda(const ConfigMap& config) {
  double omega_laser = get_laser_frequency(config);
  return (2 * MathUtils::pi / omega_laser) * PhysUtils::AtomicUnits::c;
}

inline double get_laser_period(const ConfigMap& config) {
  double omega_laser = get_laser_frequency(config);
  return (2 * MathUtils::pi / omega_laser);
}

// Number of points in the dense_frequency_spectrum=true fine frequency scan -- reads "N_omega".
// Only used by Simulation::init_simulation_parameters when dense_frequency_spectrum is true; see
// get_number_of_harmonics for the default (harmonics) mode's count, which reads "N_harmonics"
// instead -- kept as two separate keys (rather than one shared count) so switching modes doesn't
// silently reinterpret whatever count was already configured for the other mode.
inline size_t get_number_of_frequencies(const ConfigMap& config) {
  size_t N_frequencies = std::stol(get_required(config, "N_omega"));
  return N_frequencies;
}

// Number of harmonics computed in the default (dense_frequency_spectrum=false) mode -- reads
// "N_harmonics". See get_number_of_frequencies for the dense-scan mode's count ("N_omega").
inline size_t get_number_of_harmonics(const ConfigMap& config) {
  size_t N_harmonics = std::stoul(get_required(config, "N_harmonics"));
  return N_harmonics;
}

// Raw num_threads config value: 0 means "use the maximum available", resolved by the caller
// (Simulation::run_simulation) against std::thread::hardware_concurrency().
inline size_t get_num_threads(const ConfigMap& config) { return std::stoul(get_required(config, "num_threads")); }

// Returns the current user's home directory (via $HOME), so the config's output_folder can stay a
// short, machine-independent subfolder name instead of a machine-specific absolute path -- e.g. so the
// same .cfg still works after sync_to_remote.sh copies the repo to a different machine/user.
inline std::string get_home_directory() {
  const char* home = std::getenv("HOME");
  if (home == nullptr) {
    throw std::runtime_error("Cannot resolve output_folder: $HOME is not set");
  }
  return home;
}

// Creates (and returns the path to) a timestamped subfolder of "~/<output_folder>", e.g.
// "~/<output_folder>/20260709_143012", so each run's exports land in their own directory rather than
// overwriting the last run's. The config's output_folder value is always taken relative to the home
// directory (see get_home_directory above), not as an absolute path.
inline std::string make_run_output_directory(const ConfigMap& config) {
  const std::string& base_folder = get_required(config, "output_folder");

  std::time_t now_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm local_tm{};
  localtime_r(&now_time, &local_tm);

  std::ostringstream timestamp;
  timestamp << std::put_time(&local_tm, "%Y%m%d_%H%M%S");

  std::filesystem::path run_dir = std::filesystem::path(get_home_directory()) / base_folder / timestamp.str();
  std::filesystem::create_directories(run_dir);

  return run_dir.string();
}

// Copies the exact config file used for this run into its own output directory (as
// "config.cfg"), so that later inspection/plotting tools always read back the config that
// actually produced this run's data, rather than whatever config/coherent_thomson.cfg currently
// contains -- which may have since been edited to set up a different run.
inline void copy_config_to_run_directory(const std::string& config_path, const std::string& run_output_dir) {
  std::filesystem::copy_file(config_path, std::filesystem::path(run_output_dir) / "config.cfg",
                             std::filesystem::copy_options::overwrite_existing);
}

}  // namespace Core::IoUtils
