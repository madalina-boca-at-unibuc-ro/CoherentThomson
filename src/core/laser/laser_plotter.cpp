#include "../include/laser/laser_plotter.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Laser {

namespace {

// E and B extracted from F^{mu nu} (same convention as the lab-frame extraction they replace), then projected
// onto the laser's own frame (epsilon_1, epsilon_2, n), so the exported "x/y/z" components are transverse/
// transverse/longitudinal whatever laser_direction_theta/phi is. For a laser along +Oz this is the identity.
struct LaserFrameFields {
  double Ex, Ey, Ez, Bx, By, Bz;
};

LaserFrameFields to_laser_frame(const LaserField& laser, const FaradayTensor& F) {
  std::array<double, 3> E = {F[1][0], F[2][0], F[3][0]};
  std::array<double, 3> B = {F[2][3], F[3][1], F[1][2]};
  const MathUtils::RealFourVector e1 = laser.get_epsilon_1();
  const MathUtils::RealFourVector e2 = laser.get_epsilon_2();
  const MathUtils::RealFourVector n = laser.get_unity_n();
  auto project = [](const std::array<double, 3>& v, const MathUtils::RealFourVector& axis) {
    return v[0] * axis[1] + v[1] * axis[2] + v[2] * axis[3];
  };
  return {project(E, e1), project(E, e2), project(E, n), project(B, e1), project(B, e2), project(B, n)};
}

}  // namespace

void export_field_vs_phase(const LaserField& laser, const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for field export: " + filepath);
  }

  file << std::scientific << std::setprecision(6);
  file << "phi Ex Ey Ez Bx By Bz\n";  // File header

  double phi_min = laser.get_phi_min();
  double phi_max = laser.get_phi_max();
  double d_phi = laser.get_d_phi();

  // Use the physical range variables stored in your laser class instance
  double phi = phi_min;
  while (phi <= phi_max) {
    // To get fields at a specific phase, evaluate them at a position
    // where the phase equals our current target.
    // Since phi = omega*t - k*r, choosing x_mu = (c*phi/omega, 0, 0, 0) works
    // perfectly.
    double t_equivalent = phi / laser.get_omega();
    MathUtils::RealFourVector x_mu(PhysUtils::AtomicUnits::c * t_equivalent, 0.0, 0.0, 0.0);

    // Extract the Faraday Tensor
    FaradayTensor F = laser.get_faraday_tensor(x_mu);

    // Recover field components from the tensor layouts (F^{10} = Ex/c, F^{23} = -Bx, ...), in the
    // laser's own frame -- see to_laser_frame.
    auto [Ex, Ey, Ez, Bx, By, Bz] = to_laser_frame(laser, F);

    file << phi << " " << Ex << " " << Ey << " " << Ez << " " << Bx << " " << By << " " << Bz << "\n";

    phi += d_phi;
  }
}

void export_field_heatmap_z0(const LaserField& laser, double t, double x_min, double x_max, double y_min, double y_max,
                             size_t Nx, size_t Ny, double axes_scale, const std::string& axes_scale_string,
                             const std::string& filepath) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for field heatmap export: " + filepath);
  }

  double omega = laser.get_omega();
  double phi_at_t = omega * t;

  file << std::scientific << std::setprecision(6);
  file << "# axes_unit " + axes_scale_string + "\n";
  file << "# z=0 snapshot (laser frame: x/y along epsilon_1/epsilon_2, z along n) at t=" << t
       << " (phi=omega*t=" << phi_at_t << " rad = " << phi_at_t / (2 * MathUtils::pi) << " cycles)\n";
  file << "x y Ex Ey Ez Bx By Bz intensity\n";

  double dx = (Nx > 1) ? (x_max - x_min) / static_cast<double>(Nx - 1) : 0.0;
  double dy = (Ny > 1) ? (y_max - y_min) / static_cast<double>(Ny - 1) : 0.0;

  for (size_t i = 0; i < Nx; ++i) {
    double x = x_min + static_cast<double>(i) * dx;
    for (size_t j = 0; j < Ny; ++j) {
      double y = y_min + static_cast<double>(j) * dy;

      // (x, y) are coordinates in the laser's own transverse plane through the origin (its waist).
      MathUtils::RealFourVector x_mu(PhysUtils::AtomicUnits::c * t, 0.0, 0.0, 0.0);
      for (size_t mu = 1; mu < 4; ++mu) {
        x_mu[mu] = x * laser.get_epsilon_1()[mu] + y * laser.get_epsilon_2()[mu];
      }
      FaradayTensor F = laser.get_faraday_tensor(x_mu);

      // Same F^{mu nu} -> E/B extraction, in the laser's own frame, as export_field_vs_phase above.
      auto [Ex, Ey, Ez, Bx, By, Bz] = to_laser_frame(laser, F);
      double intensity = Ex * Ex + Ey * Ey + Ez * Ez;

      file << x / axes_scale << " " << y / axes_scale << " " << Ex << " " << Ey << " " << Ez << " " << Bx << " " << By
           << " " << Bz << " " << intensity << "\n";
    }
  }
}

}  // namespace Core::Laser
