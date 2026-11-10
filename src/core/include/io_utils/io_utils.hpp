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
#include <utility>
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
  } else if (to_lower(unit_name) == "omega_laser" || to_lower(unit_name) == "first_harmonic_frequency") {
    // Read directly rather than via get_laser_frequency (defined later in this file) to avoid a
    // forward reference -- same pattern the "lambda" branch above already uses.
    //
    // Both units parse to the same raw a.u. value here: the difference between "value * incident
    // laser omega" and "value * emitted (possibly Doppler-shifted) N=1 fundamental" can't be
    // resolved at config-parse time, since the fundamental depends on the dressed momentum q (built
    // from the laser and beam objects, not yet constructed when a ConfigMap alone is being parsed).
    // get_omega_range (below) is the only caller that needs this distinction; it reports which of
    // the two units was used so Simulation::init_simulation_parameters can rescale accordingly.
    return std::stod(get_required(config, "laser_frequency"));
  } else if (to_lower(unit_name) == "a.u.") {
    return 1.0;
  } else if (to_lower(unit_name) == "lambda_1") {
    // lambda_1 = 2 pi / |q_axial| needs the laser and the detector normal, which a ConfigMap alone doesn't give;
    // read_axial_length handles it for the keys that accept it. Anywhere else the fallback below would only warn
    // and use 1.0.
    throw std::invalid_argument(
        "Unit lambda_1 is accepted only by beam_cylinder_height, beam_sigma_z, beam_stack_sigma_z and "
        "beam_modulation_period");
  } else if (to_lower(unit_name) == "w0") {
    // Read directly rather than via get_laser_lg_params (defined later in this file) to avoid a
    // forward reference -- same pattern the "omega_laser" branch above uses. w0 itself is still
    // expressed in the config in an already-supported unit (typically "lambda"), not independently,
    // so this recurses through convert_unit_to_number rather than assuming atomic units directly.
    auto [w0_val, w0_unit] = split_value_and_unit(get_required(config, "laser_lg_w0"));
    return w0_val * convert_unit_to_number(w0_unit, config);
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

// Laser's propagation direction, as {theta, phi} in radians, in the simulation (lab) frame -- same
// (theta, phi) convention as get_detector_direction_angles below. Only the laser follows it: the beam
// and the detector keep their own configured lab-frame geometry.
inline std::pair<double, double> get_laser_direction_angles(const ConfigMap& config) {
  auto [theta_val, theta_unit] = split_value_and_unit(get_required(config, "laser_direction_theta"));
  double theta = theta_val * convert_unit_to_number(theta_unit, config);
  auto [phi_val, phi_unit] = split_value_and_unit(get_required(config, "laser_direction_phi"));
  double phi = phi_val * convert_unit_to_number(phi_unit, config);
  return {theta, phi};
}

// Detector's own normal direction, as {theta, phi} in radians, given in the simulation (lab) frame --
// shared by Detector::create_detector and Simulation::init_simulation_parameters so both read the same
// direction.
inline std::pair<double, double> get_detector_direction_angles(const ConfigMap& config) {
  auto [theta_val, theta_unit] = split_value_and_unit(get_required(config, "detector_direction_theta"));
  double theta = theta_val * convert_unit_to_number(theta_unit, config);
  auto [phi_val, phi_unit] = split_value_and_unit(get_required(config, "detector_direction_phi"));
  double phi = phi_val * convert_unit_to_number(phi_unit, config);
  return {theta, phi};
}

// Which frequency omega_min/omega_max's unit anchors to -- see get_omega_range below.
enum class OmegaRangeUnit { OmegaLaser, FirstHarmonicFrequency };

// omega_min/omega_max, in raw omega (not omega/c) -- only consulted when dense_frequency_spectrum
// is true (Simulation::init_simulation_parameters), to build a fine linear frequency scan instead
// of the default first-N_omega-harmonics list. Restricted to exactly two units, each with different
// downstream handling: "omega_laser" means the value is a raw multiple of the *incident* laser
// frequency and Simulation::init_simulation_parameters uses it as-is; "first_harmonic_frequency"
// means the value is a multiple of the *emitted*, possibly Doppler-shifted N=1 fundamental
// (simulation_parameters::fundamental_frequency) and needs to be rescaled by
// fundamental_frequency/laser_frequency there -- the useful choice when the beam is relativistic
// enough that the emitted frequency is far from any simple multiple of the incident one. Both units
// parse to the same numeric value here (see convert_unit_to_number); min_unit/max_unit report which
// one was actually written in the config so the caller can decide whether to rescale.
inline std::pair<double, OmegaRangeUnit> get_omega_bound(const ConfigMap& config, const std::string& key) {
  auto [val, unit] = split_value_and_unit(get_required(config, key.c_str()));
  double omega = val * convert_unit_to_number(unit, config, {"omega_laser", "first_harmonic_frequency"});
  OmegaRangeUnit unit_kind = to_lower(unit) == "first_harmonic_frequency" ? OmegaRangeUnit::FirstHarmonicFrequency
                                                                          : OmegaRangeUnit::OmegaLaser;
  return {omega, unit_kind};
}

struct OmegaRange {
  double omega_min;
  double omega_max;
  OmegaRangeUnit min_unit;
  OmegaRangeUnit max_unit;
};

inline OmegaRange get_omega_range(const ConfigMap& config) {
  auto [omega_min, min_unit] = get_omega_bound(config, "omega_min");
  auto [omega_max, max_unit] = get_omega_bound(config, "omega_max");
  return {omega_min, omega_max, min_unit, max_unit};
}

inline size_t get_laser_NT(const ConfigMap& config) { return std::stoul(get_required(config, "laser_NT")); }

inline size_t get_trajectory_NT(const ConfigMap& config) { return std::stoul(get_required(config, "trajectory_NT")); }

// Number of points in the dense_frequency_spectrum=true fine frequency scan (see get_omega_range) --
// distinct from N_harmonics (get_number_of_harmonics), which counts harmonics in the default mode.
inline size_t get_N_omega(const ConfigMap& config) { return std::stoul(get_required(config, "N_omega")); }

// The per-electron parameters a cylinder beam is generated from: cylinder geometry and Gaussian
// momentum distribution. Deliberately excludes beam_particle_count/random_seed, which control the
// generation loop rather than any individual electron. Shared by Particle::generate_cylinder_beam
// and Simulation::init_simulation_parameters (the latter only needs average_px/py/pz, for the
// initial-motion Doppler scaling of d_tau).
//
// Lengths along the beam axis (height, sigma_z, modulation.period) may be given in the unit
// lambda_1 = 2 pi / |q_axial|, which needs the laser and the detector normal. Parsing then stores them in units of
// lambda_1 (the *_in_lambda_1 flags) and Particle::resolve_beam_lengths converts them to a.u.; lengths_resolved says
// whether that has happened.
// Callers that use only the radius and the momenta (main's heatmap, init_simulation_parameters) don't need it.
struct CylinderBeamParams {
  double radius;
  bool radius_from_lg_table;  // true when beam_cylinder_radius < 0 selected the radius via lg_beam_radius_in_w0
  // H: the random cylinder's height, or the width a lattice stack's sheets fill (theory/lattice_stack-work.md, "Beam
  // distribution rework"). 0 = a flat disk.
  double height;
  // Gaussian width along the axis; 0 = sharp edges. random: beam_sigma_z, each electron displaced by N(0, sigma_z^2).
  // polar_lattice: beam_stack_sigma_z, carried by the sheets' charges (see stack_sheets). Either way the density is
  // the top-hat of width H convolved with the Gaussian.
  double sigma_z;
  bool height_in_lambda_1 = false;
  bool sigma_z_in_lambda_1 = false;
  bool lengths_resolved = true;  // false while a length above or modulation.period is still in units of lambda_1
  double lambda_1 = 0.0;         // set by resolve_beam_lengths when it is computed, 0 otherwise
  // What parsing changed or found pointless (an ignored width, coinciding sheets, a modulation that only translates):
  // printed by the beam factory and the run log, not by the parser, which runs several times per run.
  std::vector<std::string> notes;
  double center_x, center_y, center_z;
  double average_px, average_py, average_pz;
  double sigma_px, sigma_py, sigma_pz;
  // beam_distribution (optional, default random): random = independent uniform draws over the cylinder;
  // polar_lattice = deterministic equal-area rings on a flat disk, node spacing ~beam_lattice_spacing, no sampling
  // speckle (theory/future_dev/coherent_incoherent_decomposition.md section 3, step 1); beam_particle_count is then
  // ignored.
  bool polar_lattice;
  double lattice_spacing;  // requested spacing (polar_lattice only); the rings use polar_lattice_geometry's
  // Longitudinal stack (polar_lattice only, beam_stack_sheets, default 1): the same lattice repeated on M parallel
  // sheets filling the span L, cell-centred at z_m = (m - (M-1)/2) L/M, before beam_center and the beam rotation.
  // sigma_z = 0: L = H and each sheet gets 1/M of the charge. sigma_z > 0: L = H + 2 sigma_z_cut sigma_z, and each
  // sheet's charge is the smoothed density at its unmodulated z_m, normalized to a total of 1
  // (Particle::sheet_weights). Either way the total matches the flat disk.
  size_t stack_sheets;
  double sigma_z_cut = 6.0;  // beam_stack_sigma_z_cut (optional): the span's half-width beyond H/2, in sigma_z
  // Modulation (beam_modulation true, either distribution): each electron moves by its own position,
  // z -> z + (zeta / K) sin(K z + phase), K = 2 pi / period; for a random beam z includes the Gaussian wings. Charges
  // are unchanged.
  struct Modulation {
    bool active = false;
    double zeta = 0.0;    // beam_modulation_zeta: drift ratio; 0 keeps the positions
    double period = 0.0;  // beam_modulation_period: Lambda_mod, a length (a.u. once resolved)
    bool period_in_lambda_1 = false;
    double phase = 0.0;  // beam_modulation_phase (optional, default 0): phi = K z + phase
  } modulation;
};

// Rings of the polar lattice on a disk of radius R: a node at the centre and rings j = 1..J at r_j = j*d with
// round(2*pi*j) equally spaced nodes, each standing for ~d^2 of area. d is the requested spacing adjusted so that
// the outermost ring's annulus ends at R exactly, (J + 1/2) d = R. Returns {J, d}.
inline std::pair<int, double> polar_lattice_geometry(double radius, double spacing) {
  const int rings = std::max(1, static_cast<int>(std::lround(radius / spacing - 0.5)));
  return {rings, radius / (rings + 0.5)};
}

// Number of nodes on one sheet of the polar lattice with `rings` rings: the centre plus round(2*pi*j) per ring.
inline size_t polar_lattice_site_count(int rings) {
  size_t count = 1;
  for (int j = 1; j <= rings; ++j)
    count += static_cast<size_t>(std::lround(2.0 * MathUtils::pi * j));
  return count;
}

// Beam radius, in units of w0, just larger than the Laguerre-Gauss (p, l) mode's transverse extent
// (theory/table of beam electron radius in units of w_0 for different p,l.md). The mode profile depends on l only
// through |l|. Throws outside the tabulated p = 0..4, |l| = 0..6.
inline double lg_beam_radius_in_w0(int p, int l) {
  static constexpr double table[5][7] = {{1.86, 2.15, 2.37, 2.56, 2.72, 2.87, 3.01},
                                         {2.44, 2.65, 2.84, 3.00, 3.15, 3.29, 3.41},
                                         {2.87, 3.06, 3.22, 3.37, 3.50, 3.63, 3.75},
                                         {3.24, 3.40, 3.55, 3.68, 3.81, 3.93, 4.04},
                                         {3.56, 3.71, 3.84, 3.97, 4.09, 4.20, 4.31}};
  const int abs_l = std::abs(l);
  if (p < 0 || p > 4 || abs_l > 6) {
    throw std::invalid_argument("beam_cylinder_radius < 0 selects the radius from the LG table, which covers only "
                                "p = 0..4 and |l| = 0..6 (got p = " +
                                std::to_string(p) + ", l = " + std::to_string(l) + ")");
  }
  return table[p][abs_l];
}

// Picks the cylinder-beam keys (with unit conversion) out of the config into a CylinderBeamParams.
// A negative beam_cylinder_radius (its unit then ignored) means lg_beam_radius_in_w0(p, l) * w0, for a
// laguerre_gauss laser only.
// A length along the beam axis. In the unit lambda_1 it is returned in units of lambda_1, with in_lambda_1 set, for
// Particle::resolve_beam_lengths to convert; otherwise in a.u. unit_required rejects a bare number, which would be a.u.
inline double read_axial_length(const ConfigMap& config, const char* key, bool& in_lambda_1,
                                bool unit_required = false) {
  auto [value, unit] = split_value_and_unit(get_required(config, key));
  in_lambda_1 = to_lower(unit) == "lambda_1";
  if (in_lambda_1) {
    return value;
  }
  if (unit_required && unit.empty()) {
    throw std::invalid_argument(std::string(key) + " needs a unit (e.g. 1 lambda_1 or 0.5 lambda): a bare number "
                                "would be atomic units");
  }
  return value * convert_unit_to_number(unit, config);
}

inline CylinderBeamParams parse_cylinder_beam_params(const ConfigMap& config) {
  CylinderBeamParams params;

  auto [radius_val, radius_unit] = split_value_and_unit(get_required(config, "beam_cylinder_radius"));
  params.radius_from_lg_table = radius_val < 0.0;
  if (params.radius_from_lg_table) {
    if (get_required(config, "laser_type") != "laguerre_gauss") {
      throw std::invalid_argument("beam_cylinder_radius < 0 (radius from the LG table) requires "
                                  "laser_type = laguerre_gauss");
    }
    auto [p, l, w0] = get_laser_lg_params(config);
    params.radius = lg_beam_radius_in_w0(p, l) * w0;
  } else {
    params.radius = radius_val * convert_unit_to_number(radius_unit, config);
  }
  params.height = read_axial_length(config, "beam_cylinder_height", params.height_in_lambda_1);
  if (params.height < 0.0) {
    throw std::invalid_argument("beam_cylinder_height must be >= 0");
  }
  params.sigma_z = read_axial_length(config, "beam_sigma_z", params.sigma_z_in_lambda_1);
  if (params.sigma_z < 0.0) {
    throw std::invalid_argument("beam_sigma_z must be >= 0");
  }

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

  auto distribution = config.find("beam_distribution");
  const std::string distribution_name = distribution == config.end() ? "random" : distribution->second;
  if (distribution_name != "random" && distribution_name != "polar_lattice") {
    throw std::invalid_argument("beam_distribution must be random or polar_lattice (got " + distribution_name + ")");
  }
  params.polar_lattice = distribution_name == "polar_lattice";
  params.lattice_spacing = 0.0;
  if (params.polar_lattice) {
    // The lattice's Gaussian width is its own key (the type picks which keys count, below): beam_sigma_z displaces
    // electrons at random, which a lattice doesn't do.
    if (params.sigma_z != 0.0) {
      params.notes.push_back("beam_sigma_z is for random beams and is ignored by polar_lattice: use "
                             "beam_stack_sigma_z");
    }
    params.sigma_z = 0.0;
    params.sigma_z_in_lambda_1 = false;
    if (config.count("beam_stack_sigma_z") > 0) {
      params.sigma_z = read_axial_length(config, "beam_stack_sigma_z", params.sigma_z_in_lambda_1);
      if (!(params.sigma_z >= 0.0)) {
        throw std::invalid_argument("beam_stack_sigma_z must be >= 0");
      }
    }
    auto cut = config.find("beam_stack_sigma_z_cut");
    if (cut != config.end()) {
      params.sigma_z_cut = std::stod(cut->second);
      if (!(params.sigma_z_cut > 0.0)) {
        throw std::invalid_argument("beam_stack_sigma_z_cut must be > 0 (a half-width in units of sigma_z)");
      }
    }
    auto [spacing_val, spacing_unit] = split_value_and_unit(get_required(config, "beam_lattice_spacing"));
    params.lattice_spacing = spacing_val * convert_unit_to_number(spacing_unit, config);
    if (!(params.lattice_spacing > 0.0)) {
      throw std::invalid_argument("beam_lattice_spacing must be > 0");
    }
  }

  // Keys of the earlier stack layout. Unknown keys are normally ignored, but dropping one of these silently would
  // change the beam, so each names its replacement (theory/lattice_stack-work.md, "Beam distribution rework").
  const std::pair<const char*, const char*> removed_keys[] = {
      {"beam_stack_spacing", "beam_cylinder_height, the stack's width H = M d"},
      {"beam_stack_modulation", "beam_modulation true|false"},
      {"beam_modulation_sheets_per_period", "beam_stack_sheets (M) with beam_cylinder_height (H)"},
      {"beam_modulation_periods", "beam_cylinder_height (H, e.g. in lambda_1)"},
      {"beam_modulation_period_ratio", "beam_modulation_period, a length (e.g. 1 lambda_1)"}};
  for (const auto& [key, replacement] : removed_keys) {
    if (config.count(key) > 0) {
      throw std::invalid_argument(std::string(key) + " no longer exists: use " + replacement);
    }
  }

  // The type picks which keys count: a random beam ignores the lattice's (beam_lattice_spacing, beam_stack_sheets,
  // beam_stack_sigma_z, beam_stack_sigma_z_cut), as polar_lattice ignores beam_particle_count and beam_sigma_z, so
  // switching the type needs no other edit.
  params.stack_sheets = 1;
  auto stack_sheets = config.find("beam_stack_sheets");
  if (params.polar_lattice && stack_sheets != config.end()) {
    size_t parsed_length = 0;
    const long sheets = std::stol(stack_sheets->second, &parsed_length);
    if (parsed_length != stack_sheets->second.size() || sheets < 1) {
      throw std::invalid_argument("beam_stack_sheets must be an integer >= 1 (got " + stack_sheets->second + ")");
    }
    params.stack_sheets = static_cast<size_t>(sheets);
  }
  // Whether H or sigma_z is zero doesn't depend on its unit, so these cases are settled here, before lambda_1 is known.
  // With sigma_z > 0 the sheets span 2 cut sigma_z even at H = 0, so only H = 0 and sigma_z = 0 together coincide them.
  if (params.polar_lattice && params.stack_sheets == 1) {
    if (params.height > 0.0) {
      params.notes.push_back("a single sheet (beam_stack_sheets 1) has no width: beam_cylinder_height set to 0");
      params.height = 0.0;
      params.height_in_lambda_1 = false;
    }
    if (params.sigma_z > 0.0) {
      params.notes.push_back("a single sheet (beam_stack_sheets 1) can't carry a Gaussian profile: "
                             "beam_stack_sigma_z set to 0");
      params.sigma_z = 0.0;
      params.sigma_z_in_lambda_1 = false;
    }
  } else if (params.polar_lattice && params.height == 0.0 && params.sigma_z == 0.0) {
    params.notes.push_back("beam_stack_sheets " + std::to_string(params.stack_sheets) +
                           " with beam_cylinder_height 0 and beam_stack_sigma_z 0: the sheets would coincide, so the "
                           "stack is one sheet");
    params.stack_sheets = 1;
  }

  // beam_modulation (optional, default false) switches the modulation on for either distribution. Without the switch,
  // a beam_modulation_* key is an error rather than silently unused; with false, they are ignored.
  auto modulation_switch = config.find("beam_modulation");
  if (modulation_switch == config.end()) {
    for (const char* key : {"beam_modulation_zeta", "beam_modulation_period", "beam_modulation_phase"}) {
      if (config.count(key) > 0) {
        throw std::invalid_argument(std::string(key) + " needs beam_modulation true or false");
      }
    }
  } else if (modulation_switch->second != "true" && modulation_switch->second != "false") {
    throw std::invalid_argument("beam_modulation must be true or false (got " + modulation_switch->second + ")");
  }
  if (modulation_switch != config.end() && modulation_switch->second == "true") {
    CylinderBeamParams::Modulation& modulation = params.modulation;
    modulation.active = true;
    modulation.zeta = std::stod(get_required(config, "beam_modulation_zeta"));
    if (!(modulation.zeta >= 0.0)) {
      throw std::invalid_argument("beam_modulation_zeta must be >= 0");
    }
    modulation.period = read_axial_length(config, "beam_modulation_period", modulation.period_in_lambda_1, true);
    if (!(modulation.period > 0.0)) {
      throw std::invalid_argument("beam_modulation_period must be > 0");
    }
    auto phase = config.find("beam_modulation_phase");
    if (phase != config.end()) {
      auto [phase_val, phase_unit] = split_value_and_unit(phase->second);
      modulation.phase = phase_val * convert_unit_to_number(phase_unit, config);
    }
    const bool no_axial_extent =
        params.polar_lattice ? params.stack_sheets == 1 : params.height == 0.0 && params.sigma_z == 0.0;
    if (no_axial_extent) {
      params.notes.push_back("the beam has no extent along its axis, so beam_modulation only translates it: "
                             "nothing is bunched");
    }
  }

  params.lengths_resolved = !(params.height_in_lambda_1 || params.sigma_z_in_lambda_1 ||
                              (params.modulation.active && params.modulation.period_in_lambda_1));
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

// Harmonic index the default (dense_frequency_spectrum=false) mode's harmonics list starts at --
// reads "N_harmonics_min". The list then covers N_harmonics consecutive harmonic indices
// N_harmonics_min, N_harmonics_min+1, .... A value of 1 (the common case) reproduces the previous
// "first N_harmonics harmonics" behavior. Only meaningful alongside get_number_of_harmonics; unused
// in dense_frequency_spectrum=true mode, which instead scans directly between omega_min/omega_max.
inline size_t get_number_of_harmonics_min(const ConfigMap& config) {
  size_t N_harmonics_min = std::stoul(get_required(config, "N_harmonics_min"));
  return N_harmonics_min;
}

// Raw num_threads config value: 0 means "use the maximum available", resolved by the caller
// (Simulation::run_simulation) against the CPUs available to the process.
inline size_t get_num_threads(const ConfigMap& config) { return std::stoul(get_required(config, "num_threads")); }

// Returns the current user's home directory (via $HOME), so the config's output_folder can stay a
// short, machine-independent subfolder name instead of a machine-specific absolute path -- e.g. so the
// same .cfg still works on another machine/user that has the repo checked out (synced via git).
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
