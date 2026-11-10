#include "../include/particle/electron_plotter.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Core::Particle {

void plot_particle_trajectory(const std::vector<Electron>& electron_beam, const std::vector<size_t>& electron_indices,
                              const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for field export: " + filepath);
  }
  if (electron_indices.empty()) {
    throw std::runtime_error("plot_particle_trajectory: electron_indices must not be empty");
  }

  file << std::scientific << std::setprecision(6);
  // Column layout matches Electron::State: component 0 of each four-vector is
  // the time-like part, 1-3 are x/y/z.
  file << "# " << electron_indices.size() << " electron trajector" << (electron_indices.size() == 1 ? "y" : "ies")
       << "\n";
  file << "electron_id tau x0 x1 x2 x3 p0 p1 p2 p3 a0 a1 a2 a3\n";

  for (size_t electron_id : electron_indices) {
    const Electron& electron = electron_beam.at(electron_id);
    const size_t num_states = electron.get_trajectory().size();
    for (size_t i = 0; i < num_states; ++i) {
      file << electron_id << " ";
      electron.export_state(i, file);
      file << "\n";
    }
    file << "\n";  // blank line separates each electron's block
  }

  std::cout << "Successfully exported " << electron_indices.size() << " electron trajector"
            << (electron_indices.size() == 1 ? "y" : "ies") << " to " << filepath << "\n";
}

void plot_electron_beam_scatter(const std::vector<Electron>& electron_beam, const std::string& filepath,
                                double axes_scale, const std::string& axes_scale_string) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for electron beam scatter export: " + filepath);
  }

  file << std::scientific << std::setprecision(6);
  file << "# length units: " << axes_scale_string << "\n";
  // The charge weights only when they differ (a Gaussian lattice stack): equal weights don't change a normalized
  // density, and leaving the column out keeps every other beam's file as before.
  const bool weighted = std::any_of(electron_beam.begin(), electron_beam.end(), [&](const Electron& electron) {
    return electron.get_charge_weight() != electron_beam.front().get_charge_weight();
  });
  file << (weighted ? "x y z w\n" : "x y z\n");

  for (const Electron& electron : electron_beam) {
    const Core::MathUtils::RealFourVector& pos = electron.get_initial_position();
    file << pos[1] / axes_scale << " " << pos[2] / axes_scale << " " << pos[3] / axes_scale;
    if (weighted) {
      file << " " << electron.get_charge_weight();
    }
    file << "\n";
  }

  std::cout << "Successfully exported electron beam scatter positions to " << filepath << "\n";
}

}  // namespace Core::Particle