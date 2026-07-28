#pragma once

#include <string>
#include <vector>

#include "electron.hpp"

namespace Core::Particle {

void plot_particle_trajectory(const Electron& electron, const std::string& filepath);

// Exports the real (lab-frame) initial x/y/z position of every electron in the beam, to check the
// generated beam's position/spread with a 3D scatter plot.
void plot_electron_beam_scatter(const std::vector<Electron>& electron_beam, const std::string& filepath,
                                double axes_scale, const std::string& axes_scale_string);

}  // namespace Core::Particle