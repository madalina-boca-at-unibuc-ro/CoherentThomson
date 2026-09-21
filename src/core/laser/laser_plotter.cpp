#include "../include/laser/laser_plotter.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Laser {

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

    // Recover field components from the tensor layouts
    double Ex = F[1][0];  // F^{10} = Ex/c (assuming c=1 scaling)
    double Ey = F[2][0];
    double Ez = F[3][0];
    double Bx = F[2][3];  // F^{23} = -Bx
    double By = F[3][1];  // F^{31} = -By
    double Bz = F[1][2];  // F^{12} = -Bz

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
  file << "# z=0 snapshot at t=" << t << " (phi=omega*t=" << phi_at_t
       << " rad = " << phi_at_t / (2 * MathUtils::pi) << " cycles)\n";
  file << "x y Ex Ey Ez Bx By Bz intensity\n";

  double dx = (Nx > 1) ? (x_max - x_min) / static_cast<double>(Nx - 1) : 0.0;
  double dy = (Ny > 1) ? (y_max - y_min) / static_cast<double>(Ny - 1) : 0.0;

  for (size_t i = 0; i < Nx; ++i) {
    double x = x_min + static_cast<double>(i) * dx;
    for (size_t j = 0; j < Ny; ++j) {
      double y = y_min + static_cast<double>(j) * dy;

      MathUtils::RealFourVector x_mu(PhysUtils::AtomicUnits::c * t, x, y, 0.0);
      FaradayTensor F = laser.get_faraday_tensor(x_mu);

      // Same F^{mu nu} -> E/B extraction convention as export_field_vs_phase above.
      double Ex = F[1][0];
      double Ey = F[2][0];
      double Ez = F[3][0];
      double Bx = F[2][3];
      double By = F[3][1];
      double Bz = F[1][2];
      double intensity = Ex * Ex + Ey * Ey + Ez * Ez;

      file << x / axes_scale << " " << y / axes_scale << " " << Ex << " " << Ey << " " << Ez << " " << Bx << " " << By
           << " " << Bz << " " << intensity << "\n";
    }
  }
}

}  // namespace Core::Laser
