#include "../include/radiation/radiation_plotter.hpp"

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
                          const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for field export: " + filepath);
  }

  file << std::scientific << std::setprecision(6);
  file << "# coherent radiation field: one row per (frequency, screen point)\n";
  file << "# LR/SR = long_range/short_range Faraday tensor F^{mu nu}, printed as 're im' pairs\n";
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
  file << "\n";

  size_t N_omega = field.field.size();
  for (size_t i_omega = 0; i_omega < N_omega; ++i_omega) {
    size_t N_screen = field.field[i_omega].size();
    for (size_t i_screen = 0; i_screen < N_screen; ++i_screen) {
      const Simulation::Faraday& point = field.field[i_omega][i_screen];
      file << i_omega << " " << frequencies_list[i_omega] * PhysUtils::AtomicUnits::c << " "
           << i_screen;  // the frequency is stored as k = omega/c
      write_tensor(file, point.long_range);
      write_tensor(file, point.short_range);
      file << "\n";
    }
  }

  std::cout << "Successfully exported radiation field to " << filepath << "\n";
}

}  // namespace Core::Radiation
