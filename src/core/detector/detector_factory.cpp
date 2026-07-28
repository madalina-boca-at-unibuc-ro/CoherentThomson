#include "../include/detector/detector_factory.hpp"

namespace Core::Detector {

std::unique_ptr<Detector_2D>  //
create_detector(const ConfigMap& config, const Laser::LaserField& laser) {
  std::string type = IoUtils::get_required(config, "detector_type");

  if (IoUtils::get_required(config, "print_detector_type") == "true") {
    std::cout << "Creating detector of type: " << type << "\n";  // Debug output
  }

  // The detector's own normal direction, given in the canonical frame where the laser propagates along Oz
  // (independent of the laser's actual direction); rotated together with the laser below.
  auto [dir_theta, dir_phi] = IoUtils::get_detector_direction_angles(config);
  MathUtils::RealFourVector dir_vec = MathUtils::create_unit_light_like_vector<double>(dir_theta, dir_phi);
  double dir_x = dir_vec[1];
  double dir_y = dir_vec[2];
  double dir_z = dir_vec[3];

  const MathUtils::RealFourTensor laser_rotation = laser.get_rotation_matrix();

  if (type == "rectangular") {
    auto [detector_distance_val, detector_distance_unit] =
        IoUtils::split_value_and_unit(IoUtils::get_required(config, "rectangular_detector_distance"));
    detector_distance_val *= IoUtils::convert_unit_to_number(detector_distance_unit, config);
    auto [x_min_val, x_min_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "rectangular_detector_x_min"));
    x_min_val *= IoUtils::convert_unit_to_number(x_min_unit, config);
    auto [x_max_val, x_max_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "rectangular_detector_x_max"));
    x_max_val *= IoUtils::convert_unit_to_number(x_max_unit, config);
    auto [y_min_val, y_min_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "rectangular_detector_y_min"));
    y_min_val *= IoUtils::convert_unit_to_number(y_min_unit, config);
    auto [y_max_val, y_max_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "rectangular_detector_y_max"));
    y_max_val *= IoUtils::convert_unit_to_number(y_max_unit, config);

    return std::make_unique<RectangularDetector>(std::stoull(IoUtils::get_required(config, "rectangular_detector_Nx")), std::stoull(IoUtils::get_required(config, "rectangular_detector_Ny")),
                                          detector_distance_val, x_min_val, x_max_val, y_min_val, y_max_val,
                                          dir_x, dir_y, dir_z, laser_rotation);
  } else if (type == "spherical") {
    // read and include the units of the detector parameters
    auto [R_val, R_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "spherical_detector_radius"));
    R_val *= IoUtils::convert_unit_to_number(R_unit, config);
    auto [theta_min_val, theta_min_units] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "spherical_detector_theta_min"));
    theta_min_val *= IoUtils::convert_unit_to_number(theta_min_units, config);
    auto [theta_max_val, theta_max_units] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "spherical_detector_theta_max"));
    theta_max_val *= IoUtils::convert_unit_to_number(theta_max_units, config);
    auto [phi_min_val, phi_min_units] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "spherical_detector_phi_min"));
    phi_min_val *= IoUtils::convert_unit_to_number(phi_min_units, config);
    auto [phi_max_val, phi_max_units] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "spherical_detector_phi_max"));
    phi_max_val *= IoUtils::convert_unit_to_number(phi_max_units, config);

    return std::make_unique<SphericalDetector>(
        std::stoull(IoUtils::get_required(config, "spherical_detector_N_theta")), std::stoull(IoUtils::get_required(config, "spherical_detector_N_phi")),
        R_val,  // already converted to atomic units above via convert_unit_to_number(R_unit, config)
        theta_min_val, theta_max_val, phi_min_val, phi_max_val, dir_x, dir_y, dir_z, laser_rotation);
  } else if (type == "circular") {
    auto [detector_distance_val, detector_distance_unit] =
        IoUtils::split_value_and_unit(IoUtils::get_required(config, "circular_detector_distance"));
    detector_distance_val *= IoUtils::convert_unit_to_number(detector_distance_unit, config);
    auto [R_min_val, R_min_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "circular_detector_R_min"));
    R_min_val *= IoUtils::convert_unit_to_number(R_min_unit, config);
    auto [R_max_val, R_max_unit] = IoUtils::split_value_and_unit(IoUtils::get_required(config, "circular_detector_R_max"));
    R_max_val *= IoUtils::convert_unit_to_number(R_max_unit, config);

    return std::make_unique<CircularDetector>(
        std::stoull(IoUtils::get_required(config, "circular_detector_N_R")), std::stoull(IoUtils::get_required(config, "circular_detector_N_phi")),
        detector_distance_val, R_min_val, R_max_val, dir_x, dir_y, dir_z, laser_rotation);
  } else {
    throw std::runtime_error("Unknown screen type: " + type);
  }
}
}  // namespace Core::Detector