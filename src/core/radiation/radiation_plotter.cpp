#include "../include/radiation/radiation_plotter.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

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
          "is identically zero for radiation_formula=\"direct\" (see theory/"
          "FT_Faraday_tensor-direct_and_simplified_forms.md's \"Form 2's boundary term F_b\" section)\n";
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

void export_incident_field_fourier(const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                                   double fundamental_frequency, const std::string& filepath) {
  std::string detector_type = detector.get_type_name();
  if (detector_type != "RectangularDetector" && detector_type != "CircularDetector") {
    throw std::runtime_error("export_incident_field_fourier: detector_type '" + detector_type +
                             "' has no well-defined flat local (x, y) plane -- only RectangularDetector and "
                             "CircularDetector are supported (matching the angular-momentum theory docs' own "
                             "flat-screen restriction).");
  }

  // Middle of the pulse's flat-top plateau (envelope = 0 there, i.e. peak amplitude) -- see this
  // function's own doc comment (radiation_plotter.hpp) for why the exact instant chosen doesn't matter
  // to any of the downstream bilinear angular-momentum formulas.
  double phi_mid = 0.5 * (laser.get_phi_min() + laser.get_phi_max());
  double ct0 = phi_mid * PhysUtils::AtomicUnits::c / laser.get_omega();
  Core::MathUtils::RealFourVector eps1 = laser.get_epsilon_1();
  Core::MathUtils::RealFourVector eps2 = laser.get_epsilon_2();

  size_t N1 = detector.get_cols();
  size_t N2 = detector.get_rows();

  Simulation::RadiationField incident_field;
  incident_field.field.assign(1, std::vector<Simulation::Faraday>(detector.get_total_points()));

  for (size_t i = 0; i < N1; ++i) {
    for (size_t j = 0; j < N2; ++j) {
      double x_loc, y_loc;
      if (detector_type == "CircularDetector") {
        // get_row_coordinate/get_col_coordinate return (r, phi) for a circular detector, not (x, y)
        // directly -- see Detector::CircularDetector's own constructor for the identical formula.
        double r = detector.get_row_coordinate(i);
        double phi_angle = detector.get_col_coordinate(j);
        x_loc = r * std::cos(phi_angle);
        y_loc = r * std::sin(phi_angle);
      } else {
        x_loc = detector.get_row_coordinate(i);
        y_loc = detector.get_col_coordinate(j);
      }

      // Canonical-frame position (t=ct0/c, x_loc, y_loc, z_loc=0), built directly from the laser's own
      // transverse basis vectors -- bypassing the detector's actual to_lab_frame rotation/distance
      // entirely, so this screen sits at the laser's own waist regardless of the real detector geometry.
      Core::MathUtils::RealFourVector x_mu;
      x_mu[0] = ct0;
      for (size_t mu = 1; mu < 4; ++mu) {
        x_mu[mu] = x_loc * eps1[mu] + y_loc * eps2[mu];
      }

      size_t i_screen = i * N2 + j;
      incident_field.field[0][i_screen].long_range = laser.get_complex_faraday_tensor(x_mu);
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
