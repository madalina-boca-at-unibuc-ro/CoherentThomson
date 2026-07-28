#pragma once

#include <string>

#include "detector.hpp"

namespace Core::Detector {

void plot_detector(const Detector_2D& detector, const std::string& type_of_plot, const std::string& output_dir,
                   const double& axes_scale, const std::string& axes_scale_string);

// Exports the real (lab-frame) x/y/z position of every detector point, regardless of detector type, so
// the detector's actual position and orientation can be checked with a 3D scatter plot.
void plot_detector_scatter(const Detector_2D& detector, const std::string& output_dir, const double& axes_scale,
                           const std::string& axes_scale_string);

}  // namespace Core::Detector