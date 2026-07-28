#pragma once
#include <cstddef>
#include <string>

#include "laser_field.hpp"

namespace Core::Laser {
void export_field_vs_phase(const LaserField& laser, const std::string& filepath);

// Exports a heat map of the field intensity (Ex^2+Ey^2+Ez^2) over an x/y grid, in the laser's own
// canonical frame (laser along Oz) at z=0 and a fixed time t -- i.e. independent of the laser's actual
// configured laser_nx/ny/nz direction. x_min/x_max/y_min/y_max are given in the same atomic units as
// the rest of the physics (already unit-converted by the caller); axes_scale/axes_scale_string follow
// the same "divide before writing, label the unit in a header comment" convention as
// Detector::plot_detector_scatter.
void export_field_heatmap_z0(const LaserField& laser, double t, double x_min, double x_max, double y_min, double y_max,
                             size_t Nx, size_t Ny, double axes_scale, const std::string& axes_scale_string,
                             const std::string& filepath);
}  // namespace Core::Laser
