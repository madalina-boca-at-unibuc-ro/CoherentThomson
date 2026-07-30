#pragma once

#include <string>
#include <vector>

#include "electron.hpp"

namespace Core::Particle {

// Exports the trajectories of the given electrons (identified by their index within
// electron_beam) into one file, each row tagged with its electron_id column so a single
// plotting pass can separate them back out.
void plot_particle_trajectory(const std::vector<Electron>& electron_beam,
                              const std::vector<size_t>& electron_indices, const std::string& filepath);

// Exports the real (lab-frame) initial x/y/z position of every electron in the beam, to check the
// generated beam's position/spread with a 3D scatter plot.
void plot_electron_beam_scatter(const std::vector<Electron>& electron_beam, const std::string& filepath,
                                double axes_scale, const std::string& axes_scale_string);

}  // namespace Core::Particle