#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "../io_utils/config_map.hpp"
#include "../io_utils/io_utils.hpp"
#include "../laser/laser_field.hpp"
#include "../math_utils/math_constants.hpp"
#include "../phys_utils/phys_utils.hpp"
#include "detector.hpp"

namespace Core::Detector {

// laser's spatial rotation matrix is shared with the detector so it's built in its own canonical-frame
// direction (detector_direction_theta/phi) and then carried along with the laser's actual orientation.
std::unique_ptr<Detector_2D> create_detector(const ConfigMap& config, const Laser::LaserField& laser);
}  // namespace Core::Detector