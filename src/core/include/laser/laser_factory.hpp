#pragma once
#include <memory>
#include <stdexcept>
#include <string>

#include "../io_utils/config_map.hpp"
#include "../io_utils/io_utils.hpp"
#include "../math_utils/math_utils.hpp"
#include "laser_field.hpp"

namespace Core::Laser {

std::unique_ptr<LaserField> create_laser(const ConfigMap& config);
}