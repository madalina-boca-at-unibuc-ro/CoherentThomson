#include "../include/laser/laser_factory.hpp"

namespace Core::Laser {

std::unique_ptr<LaserField> create_laser(const ConfigMap& config) {
  try {
    // 1. Parse scalar configuration parameters
    double omega = IoUtils::get_laser_frequency(config);
    double a0 = IoUtils::get_laser_a0(config);
    double flat_duration_val = IoUtils::get_laser_flat_duration(config);
    double wing_sigma_val = IoUtils::get_laser_wing_sigma(config);
    double delay_val = IoUtils::get_laser_delay(config);
    double wing_cutoff_sigmas = IoUtils::get_laser_wing_sigma_cutoff(config);
    // TO DO: the delay here should be transferred to the initial time in simulation : to be written in CLAUDE.md

    // 2. Parse the 3D propagation direction vector components
    Core::MathUtils::RealFourVector unit_n_in = IoUtils::get_laser_direction(config);  // Temporarily store as a
                                                                                       // 4-vector for normalization

    auto [zeta_1_in, zeta_2_in] = IoUtils::get_laser_zeta(config);

    size_t NT_in = IoUtils::get_laser_NT(config);

    // 3. Dispatch on the laser's wave type
    std::string laser_type = IoUtils::get_required(config, "laser_type");
    if (laser_type == "plane_wave") {
      return std::make_unique<PlaneWaveLaser>(omega, a0, flat_duration_val, wing_sigma_val, delay_val,
                                              wing_cutoff_sigmas, zeta_1_in, zeta_2_in, NT_in, unit_n_in);
    } else if (laser_type == "laguerre_gauss") {
      auto [p_val, l_val, w0_val] = IoUtils::get_laser_lg_params(config);
      return std::make_unique<LaguerreGaussLaser>(omega, a0, flat_duration_val, wing_sigma_val, delay_val,
                                                  wing_cutoff_sigmas, zeta_1_in, zeta_2_in, NT_in, unit_n_in, p_val,
                                                  l_val, w0_val);
    } else {
      throw std::runtime_error("LaserFactory Error: Unknown laser_type \"" + laser_type + "\"");
    }
  } catch (const std::out_of_range& e) {
    throw std::runtime_error("LaserFactory Error: Missing required laser configuration key! " + std::string(e.what()));
  } catch (const std::invalid_argument& e) {
    throw std::runtime_error("LaserFactory Error: Invalid numerical format in "
                             "laser configuration values!");
  }
}
}  // namespace Core::Laser