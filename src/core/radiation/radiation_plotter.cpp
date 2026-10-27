#include "../include/radiation/radiation_plotter.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../include/detector/detector.hpp"
#include "../include/io_utils/io_utils.hpp"
#include "../include/laser/laser_field.hpp"
#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Radiation {

namespace {
// Writes all 16 F^{mu nu} elements of a complex FourTensor as "re im" pairs, mu-major, nu-minor.
// Leading (not trailing) separators, so callers never leave a stray trailing space before the
// newline -- a trailing space would tokenize into a phantom extra field and misalign every
// column against the header when read back with a whitespace-delimited parser (e.g. pandas).
void write_tensor(std::ofstream& file, const MathUtils::ComplexFourTensor& tensor) {
  auto components = tensor.uu();
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = 0; nu < 4; ++nu) {
      file << " " << components[mu][nu].real() << " " << components[mu][nu].imag();
    }
  }
}
}  // namespace

void plot_radiation_field(const Simulation::RadiationField& field, const std::vector<double>& frequencies_list,
                          double fundamental_frequency, const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for field export: " + filepath);
  }

  file << std::scientific << std::setprecision(6);
  file << "# coherent radiation field: one row per (frequency, screen point)\n";
  file << "# LR/SR/BR = long_range/short_range/boundary Faraday tensor F^{mu nu}, printed as 're im' pairs -- BR "
          "is identically zero for radiation_formula=\"direct\"/\"long_distance_direct\", SR for both "
          "long_distance_* formulas; only LR+SR+BR is physical\n";
  file << "# omega is in units of the fundamental (non_linear_Thomson_formula(k1, p, n2, 1)); 1.0 = fundamental, "
          "3.0 = third harmonic\n";
  file << "i_omega omega i_screen";
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = 0; nu < 4; ++nu) {
      file << " LR_F" << mu << nu << "_re LR_F" << mu << nu << "_im";
    }
  }
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = 0; nu < 4; ++nu) {
      file << " SR_F" << mu << nu << "_re SR_F" << mu << nu << "_im";
    }
  }
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t nu = 0; nu < 4; ++nu) {
      file << " BR_F" << mu << nu << "_re BR_F" << mu << nu << "_im";
    }
  }
  file << "\n";

  size_t N_omega = field.field.size();
  for (size_t i_omega = 0; i_omega < N_omega; ++i_omega) {
    size_t N_screen = field.field[i_omega].size();
    for (size_t i_screen = 0; i_screen < N_screen; ++i_screen) {
      const Simulation::Faraday& point = field.field[i_omega][i_screen];
      file << i_omega << " " << frequencies_list[i_omega] / fundamental_frequency << " "
           << i_screen;  // omega in units of the fundamental (both frequencies_list and fundamental_frequency are k =
                         // omega/c, so the ratio is dimensionless)
      write_tensor(file, point.long_range);
      write_tensor(file, point.short_range);
      write_tensor(file, point.boundary);
      file << "\n";
    }
  }

  std::cout << "Successfully exported radiation field to " << filepath << "\n";
}

namespace {

// F'^{ab} = Lambda^a_k Lambda^b_l F^{kl}, with Lambda's spatial rows the laser-frame axes (epsilon_1, epsilon_2, n)
// in lab components and Lambda^0_0 = 1: the lab tensor's components in the laser's own frame. Identity for a laser
// along +Oz (epsilon_1 = x, epsilon_2 = y, n = z).
MathUtils::ComplexFourTensor to_laser_frame(const MathUtils::ComplexFourTensor& F, const MathUtils::RealFourVector& e1,
                                            const MathUtils::RealFourVector& e2, const MathUtils::RealFourVector& n) {
  std::array<std::array<double, 4>, 4> lambda{};
  lambda[0][0] = 1.0;
  for (size_t k = 1; k < 4; ++k) {
    lambda[1][k] = e1[k];
    lambda[2][k] = e2[k];
    lambda[3][k] = n[k];
  }
  MathUtils::ComplexFourTensor result;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      MathUtils::Complex sum = 0.0;
      for (size_t k = 0; k < 4; ++k) {
        for (size_t l = 0; l < 4; ++l) {
          sum += lambda[a][k] * lambda[b][l] * F[k][l];
        }
      }
      result[a][b] = sum;
    }
  }
  return result;
}

}  // namespace

IncidentScreen incident_screen(const Laser::LaserField& laser, const Detector::Detector_2D& detector) {
  IncidentScreen screen{0.0, 0.0};
  if (const auto* circular = dynamic_cast<const Detector::CircularDetector*>(&detector)) {
    screen.distance = circular->get_distance();
  } else if (const auto* rectangular = dynamic_cast<const Detector::RectangularDetector*>(&detector)) {
    screen.distance = rectangular->get_distance();
  } else {
    throw std::runtime_error("incident screen: detector_type '" + detector.get_type_name() +
                             "' has no well-defined flat local (x, y) plane -- only RectangularDetector and "
                             "CircularDetector are supported.");
  }
  // The plane is downstream of the waist along the laser whichever way the detector faces.
  screen.distance = std::abs(screen.distance);
  if (const auto* lg = dynamic_cast<const Laser::LaguerreGaussLaser*>(&laser)) {
    double z_ratio = screen.distance / lg->get_lg_rayleigh_range();
    double w_Z = lg->get_lg_w0() * std::sqrt(1.0 + z_ratio * z_ratio);
    screen.half_width = 2.0 * IoUtils::lg_beam_radius_in_w0(lg->get_lg_p(), lg->get_lg_l()) * w_Z;
  }
  return screen;
}

void write_incident_config(const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                           const std::string& run_output_dir) {
  IncidentScreen screen = incident_screen(laser, detector);
  bool circular = detector.get_type_name() == "CircularDetector";
  std::string prefix = circular ? "circular_detector_" : "rectangular_detector_";
  std::filesystem::path run_dir(run_output_dir);
  std::ofstream file(run_dir / "incident_config.cfg");
  if (!file) {
    throw std::runtime_error("Could not open " + (run_dir / "incident_config.cfg").string() + " for writing");
  }

  file << "# Written by Core::Radiation::write_incident_config: the incident_field.dat screen, for the Python\n"
       << "# scripts' --incident (they read a key's first occurrence, so these lines override the copy of\n"
       << "# config.cfg below). Laser frame, plane z = distance downstream of the waist. Not a solver config.\n";
  file << std::setprecision(17);
  file << "detector_direction_theta " << 0.0 << " pi\n"
       << "detector_direction_phi " << 0.0 << " pi\n"
       << prefix << "distance " << screen.distance << " a.u.\n";
  if (screen.half_width > 0.0) {
    if (circular) {
      file << "circular_detector_R_min " << 0.0 << " a.u.\n"
           << "circular_detector_R_max " << screen.half_width << " a.u.\n";
    } else {
      for (const char* axis : {"x", "y"}) {
        file << "rectangular_detector_" << axis << "_min " << -screen.half_width << " a.u.\n"
             << "rectangular_detector_" << axis << "_max " << screen.half_width << " a.u.\n";
      }
    }
  }
  file << "# ---------------- copy of config.cfg ----------------\n";
  std::ifstream original(run_dir / "config.cfg");
  file << original.rdbuf();
}

void export_incident_field_fourier(const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                                   double fundamental_frequency, const std::string& filepath) {
  std::string detector_type = detector.get_type_name();
  IncidentScreen screen = incident_screen(laser, detector);  // throws for detectors without a flat plane

  // Middle of the pulse's flat-top plateau (envelope = 0 there, i.e. peak amplitude) -- see this
  // function's own doc comment (radiation_plotter.hpp) for why the exact instant chosen doesn't matter.
  // ct is shifted by the plane's distance so that phi = k (ct - z_loc) is the same as at the waist.
  double phi_mid = 0.5 * (laser.get_phi_min() + laser.get_phi_max());
  double ct0 = phi_mid * PhysUtils::AtomicUnits::c / laser.get_omega() + screen.distance;
  Core::MathUtils::RealFourVector eps1 = laser.get_epsilon_1();
  Core::MathUtils::RealFourVector eps2 = laser.get_epsilon_2();
  Core::MathUtils::RealFourVector n = laser.get_unity_n();

  size_t N1 = detector.get_cols();
  size_t N2 = detector.get_rows();
  // Own grid (half_width > 0): the detector's construction (CircularDetector's equal-area r_i with R_min = 0,
  // RectangularDetector's linspace) over the beam-sized window, with the detector's point counts.
  double R = screen.half_width;
  auto own_coordinate = [&](size_t i, size_t count) {
    if (count < 2) return 0.0;
    double fraction = static_cast<double>(i) / static_cast<double>(count - 1);
    return detector_type == "CircularDetector" ? R * std::sqrt(fraction) : -R + 2.0 * R * fraction;
  };

  Simulation::RadiationField incident_field;
  incident_field.field.assign(1, std::vector<Simulation::Faraday>(detector.get_total_points()));

  for (size_t i = 0; i < N1; ++i) {
    for (size_t j = 0; j < N2; ++j) {
      double x_loc, y_loc;
      if (detector_type == "CircularDetector") {
        // get_row_coordinate/get_col_coordinate return (r, phi) for a circular detector, not (x, y)
        // directly -- see Detector::CircularDetector's own constructor for the identical formula.
        double r = R > 0.0 ? own_coordinate(i, N1) : detector.get_row_coordinate(i);
        double phi_angle = detector.get_col_coordinate(j);
        x_loc = r * std::cos(phi_angle);
        y_loc = r * std::sin(phi_angle);
      } else {
        x_loc = R > 0.0 ? own_coordinate(i, N1) : detector.get_row_coordinate(i);
        y_loc = R > 0.0 ? own_coordinate(j, N2) : detector.get_col_coordinate(j);
      }

      // Position (t=ct0/c, x_loc, y_loc, z_loc=Z) in the laser's own frame, built directly from its basis
      // vectors -- bypassing the detector's actual to_lab_frame rotation entirely. The field components are
      // written in the same laser frame (see to_laser_frame), so the file shows the laser as if it propagated
      // along +Oz, whatever laser_direction_theta/phi is.
      Core::MathUtils::RealFourVector x_mu;
      x_mu[0] = ct0;
      for (size_t mu = 1; mu < 4; ++mu) {
        x_mu[mu] = x_loc * eps1[mu] + y_loc * eps2[mu] + screen.distance * n[mu];
      }

      size_t i_screen = i * N2 + j;
      incident_field.field[0][i_screen].long_range = to_laser_frame(laser.get_complex_faraday_tensor(x_mu), eps1,
                                                                    eps2, laser.get_unity_n());
      // short_range and boundary stay default-constructed (identically zero): the incident field has
      // no such split, it's simply "the whole field" -- see this function's own doc comment.
    }
  }

  // fundamental_frequency (like every entry of sim_par.frequencies) is in PhysUtils::non_linear_Thomson_formula's
  // omega/c convention, not raw angular frequency -- frequencies_list here must match that convention
  // for plot_radiation_field's "frequencies_list[i]/fundamental_frequency" division to come out
  // dimensionless (omega/omega_1) instead of picking up a spurious factor of c.
  double k_wave = laser.get_omega() / PhysUtils::AtomicUnits::c;

  // Reuses plot_radiation_field verbatim for the actual file writing, guaranteeing an identical format
  // to radiation_field.dat -- its own stdout line says "radiation field", so print a clearer one after.
  plot_radiation_field(incident_field, {k_wave}, fundamental_frequency, filepath);
  std::cout << "Successfully exported incident field to " << filepath << "\n";
}

}  // namespace Core::Radiation
