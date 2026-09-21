#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "../io_utils/config_map.hpp"
#include "../io_utils/io_utils.hpp"
#include "../math_utils/math_constants.hpp"
#include "../phys_utils/phys_utils.hpp"
#include "detector.hpp"

namespace Core::Detector {

// Built in its own direction (detector_direction_theta/phi), within the simulation frame (the laser
// always propagates along Oz).
std::unique_ptr<Detector_2D> create_detector(const ConfigMap& config);
}  // namespace Core::Detector