#include "../include/detector/detector.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Core::Detector {

// =========================================================================
// BASE CLASS DETECTOR_2D DEFINITIONS (This fixes the missing vtable!)
// =========================================================================

Detector_2D::Detector_2D(size_t n1, size_t n2, double dir_x, double dir_y, double dir_z,
                         const MathUtils::RealFourTensor& laser_rotation)
    : N1(n1),
      N2(n2),
      local_rotation(MathUtils::rotation_matrix_from_direction(dir_x, dir_y, dir_z)),
      lab_rotation(laser_rotation) {
  points.reserve(n1 * n2);
  MathUtils::RealFourVector normal_dir_4 = to_lab_frame(0.0, 0.0, 1.0);
  normal_dir = {normal_dir_4[1], normal_dir_4[2], normal_dir_4[3]};
}

Detector_2D::~Detector_2D() {}  // Explicitly defines the virtual destructor anchor point

MathUtils::RealFourVector Detector_2D::to_lab_frame(double x_local, double y_local, double z_local) const {
  std::array<double, 3> p_canonical = MathUtils::rotate3d(local_rotation, {x_local, y_local, z_local});
  MathUtils::RealFourVector v_canonical(0.0, p_canonical[0], p_canonical[1], p_canonical[2]);
  return MathUtils::contract(lab_rotation, v_canonical);
}

void Detector_2D::print_info() const {
  std::cout << "========================================\n"
            << "       DETECTOR SCREEN DIAGNOSTICS      \n"
            << "========================================\n"
            << "Grid Resolution : " << N1 << " x " << N2 << "\n"
            << "Total Points    : " << get_total_points() << "\n";
}

// =========================================================================
// SPHERICAL DETECTOR IMPLEMENTATION
// =========================================================================

SphericalDetector::SphericalDetector(size_t N_theta, size_t N_phi, double R, double theta_cone_min,
                                     double theta_cone_max, double phi_cone_min, double phi_cone_max, double dir_x,
                                     double dir_y, double dir_z, const MathUtils::RealFourTensor& laser_rotation)
    : Detector_2D(N_theta, N_phi, dir_x, dir_y, dir_z, laser_rotation),
      radius(R),
      phi_cone_min(phi_cone_min),
      phi_cone_max(phi_cone_max) {
  cos_theta_cone_min = std::cos(theta_cone_min);
  cos_theta_cone_max = std::cos(theta_cone_max);
  d_cos_theta_cone = (N_theta > 1) ? (cos_theta_cone_max - cos_theta_cone_min) / static_cast<double>(N_theta - 1) : 0.0;
  d_phi_cone = (N_phi > 1) ? (phi_cone_max - phi_cone_min) / static_cast<double>(N_phi - 1) : 0.0;

  for (size_t i = 0; i < N_theta; ++i) {
    double cos_theta = cos_theta_cone_min + static_cast<double>(i) * d_cos_theta_cone;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    for (size_t j = 0; j < N_phi; ++j) {
      double phi = phi_cone_min + static_cast<double>(j) * d_phi_cone;

      // local position on a cone along oz
      double x_local = radius * sin_theta * std::cos(phi);
      double y_local = radius * sin_theta * std::sin(phi);
      double z_local = radius * cos_theta;

      points.push_back(to_lab_frame(x_local, y_local, z_local));
    }
  }
}

double SphericalDetector::get_row_coordinate(size_t i) const {
  return cos_theta_cone_min + static_cast<double>(i) * d_cos_theta_cone;
}

double SphericalDetector::get_col_coordinate(size_t j) const {
  return phi_cone_min + static_cast<double>(j) * d_phi_cone;
}

void SphericalDetector::print_info() const {
  // 1. Call the parent function first to print the shared grid stats
  Detector_2D::print_info();

  // 2. Print only what makes a spherical detector unique
  std::cout << "Screen Geometry : Spherical Mesh\n"
            << "Radius          : " << radius << "\n"
            << "========================================\n";
}

// =========================================================================
// RECTANGULAR DETECTOR IMPLEMENTATION
// =========================================================================

RectangularDetector::RectangularDetector(size_t Nx, size_t Ny, double D, double x1, double x2, double y1, double y2,
                                         double dir_x, double dir_y, double dir_z,
                                         const MathUtils::RealFourTensor& laser_rotation)
    : Detector_2D(Nx, Ny, dir_x, dir_y, dir_z, laser_rotation),
      distance(D),
      x_min(x1),
      x_max(x2),
      y_min(y1),
      y_max(y2) {
  dx = (Nx > 1) ? (x_max - x_min) / static_cast<double>(Nx - 1) : 0.0;
  dy = (Ny > 1) ? (y_max - y_min) / static_cast<double>(Ny - 1) : 0.0;

  // Fill the 1D array in the parent
  for (size_t i = 0; i < N1; ++i) {
    for (size_t j = 0; j < N2; ++j) {
      double x_local = x_min + static_cast<double>(i) * dx;
      double y_local = y_min + static_cast<double>(j) * dy;
      points.push_back(to_lab_frame(x_local, y_local, distance));
    }
  }
}

double RectangularDetector::get_row_coordinate(size_t i) const { return x_min + static_cast<double>(i) * dx; }

double RectangularDetector::get_col_coordinate(size_t j) const { return y_min + static_cast<double>(j) * dy; }

void RectangularDetector::print_info() const {
  // 1. Call the parent function first to print the shared grid stats
  Detector_2D::print_info();

  // 2. Print only what makes a rectangular detector unique
  std::cout << "Screen Geometry : Rotatable Rectangular Plane\n"
            << "Distance (D)    : " << distance << "\n"
            << "========================================\n";
}

// =========================================================================
// CIRCULAR DETECTOR IMPLEMENTATION
// =========================================================================

CircularDetector::CircularDetector(size_t N_R, size_t N_phi, double D, double r_min, double r_max, double dir_x,
                                   double dir_y, double dir_z, const MathUtils::RealFourTensor& laser_rotation)
    : Detector_2D(N_R, N_phi, dir_x, dir_y, dir_z, laser_rotation), distance(D), R_min(r_min), R_max(r_max) {
  // Equal-area radial spacing (r_i^2 linear in i), not equal-distance -- see the class comment in detector.hpp.
  dR_sq = (N_R > 1) ? (R_max * R_max - R_min * R_min) / static_cast<double>(N_R - 1) : 0.0;
  d_phi = (N_phi > 1) ? 2.0 * MathUtils::pi / static_cast<double>(N_phi - 1) : 0.0;

  for (size_t i = 0; i < N1; ++i) {
    double r = std::sqrt(R_min * R_min + static_cast<double>(i) * dR_sq);

    for (size_t j = 0; j < N2; ++j) {
      double phi = static_cast<double>(j) * d_phi;

      // local position on a flat annulus at z_local = distance
      double x_local = r * std::cos(phi);
      double y_local = r * std::sin(phi);

      points.push_back(to_lab_frame(x_local, y_local, distance));
    }
  }
}

double CircularDetector::get_row_coordinate(size_t i) const {
  return std::sqrt(R_min * R_min + static_cast<double>(i) * dR_sq);
}

double CircularDetector::get_col_coordinate(size_t j) const { return static_cast<double>(j) * d_phi; }

void CircularDetector::print_info() const {
  // 1. Call the parent function first to print the shared grid stats
  Detector_2D::print_info();

  // 2. Print only what makes a circular detector unique
  std::cout << "Screen Geometry : Circular Annulus\n"
            << "R_min           : " << R_min << "\n"
            << "R_max           : " << R_max << "\n"
            << "Distance (D)    : " << distance << "\n"
            << "========================================\n";
}

void Detector_2D::export_coordinates(const std::string& filepath, double scale, std::string scale_units) const {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file for exporting screen coordinates: " + filepath);
  }

  // Apply clean scientific formatting
  file << std::scientific << std::setprecision(6);

  file << "# length units: " << scale_units << "\n";
  file << "x y z\n";

  // Loop over internal grid dimensions directly
  for (const auto& pos : points) {
    // Write spatial coordinates (indices 1, 2, 3 of our 4-vector)
    file << pos[1] / scale << " " << pos[2] / scale << " " << pos[3] / scale << "\n";
  }

  // File automatically closes when the std::ofstream object goes out of scope
  // here
}

std::pair<double, double> SphericalDetector::get_stereographic_projection(size_t index) const {
  auto [i, j] = get_grid_indices(index);
  double cos_theta = get_row_coordinate(i);
  double phi = get_col_coordinate(j);
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double r_proj = sin_theta / (1.0 + cos_theta);
  double rho = radius * r_proj;
  double x_proj = rho * std::cos(phi);
  double y_proj = rho * std::sin(phi);
  return {x_proj, y_proj};
}

}  // namespace Core::Detector