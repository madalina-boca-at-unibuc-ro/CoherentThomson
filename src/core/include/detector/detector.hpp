#pragma once
#include <array>
#include <utility>
#include <vector>

#include "../include/math_utils/math_utils.hpp"

namespace Core::Detector {

class Detector_2D {
protected:
  size_t N1, N2;  // Grid dimensions (Nx, Ny) or (N_theta, N_phi)
  std::vector<MathUtils::RealFourVector> points;

  // The detector's own local rotation, orthogonal to its canonical-frame direction (dir_x/y/z), and the laser's
  // full 4x4 rotation tensor that carries it into the lab frame -- kept as two physically distinct steps rather
  // than pre-combined into one matrix.
  MathUtils::RotationMatrix3x3 local_rotation;
  MathUtils::RealFourTensor lab_rotation;
  std::array<double, 3> normal_dir;  // detector's own normal direction, in the lab frame

  // Rotates a point given in the detector's local canonical frame into the lab frame: first the local 3x3
  // rotation orthogonal to the detector's own direction, then the laser's full 4x4 rotation tensor.
  MathUtils::RealFourVector to_lab_frame(double x_local, double y_local, double z_local) const;

public:
  // dir_x/y/z is the detector's own normal direction in the canonical frame (as if the laser pointed along Oz).
  Detector_2D(size_t n1, size_t n2, double dir_x, double dir_y, double dir_z,
              const MathUtils::RealFourTensor& laser_rotation);
  virtual ~Detector_2D();

  size_t get_cols() const { return N1; }
  size_t get_rows() const { return N2; }
  size_t get_total_points() const { return N1 * N2; }
  MathUtils::RealFourVector get_point(size_t ix) const { return points.at(ix); }
  const std::array<double, 3>& get_normal_direction() const { return normal_dir; }

  std::pair<size_t, size_t> get_grid_indices(size_t index) const { return {index / N2, index % N2}; }

  // Pure virtual function: This handles the one-to-one mapping to 3D space
  virtual std::string get_type_name() const = 0;
  virtual double get_row_coordinate(size_t i) const = 0;  // For diagnostic output
  virtual double get_col_coordinate(size_t j) const = 0;  // For diagnostic output
  virtual void print_info() const;
  void export_coordinates(const std::string& filepath, double scale, std::string scale_units) const;
};

class SphericalDetector : public Detector_2D {
private:
  double radius;
  double cos_theta_cone_min;
  double cos_theta_cone_max;
  double phi_cone_min;
  double phi_cone_max;
  double d_cos_theta_cone, d_phi_cone;

public:
  SphericalDetector(size_t N_theta, size_t N_phi, double R, double theta_cone_min, double theta_cone_max,
                    double phi_cone_min, double phi_cone_max, double dir_x, double dir_y, double dir_z,
                    const MathUtils::RealFourTensor& laser_rotation);
  std::string get_type_name() const override { return "SphericalDetector"; }
  double get_row_coordinate(size_t i) const override;
  double get_col_coordinate(size_t j) const override;
  void print_info() const override;
  std::pair<double, double> get_stereographic_projection(size_t index) const;
};

class RectangularDetector : public Detector_2D {
private:
  double distance;
  double x_min, x_max;
  double y_min, y_max;
  double dx, dy;

public:
  RectangularDetector(size_t Nx, size_t Ny, double D, double x1, double x2, double y1, double y2, double dir_x,
                      double dir_y, double dir_z, const MathUtils::RealFourTensor& laser_rotation);

  std::string get_type_name() const override { return "RectangularDetector"; }
  double get_row_coordinate(size_t i) const override;
  double get_col_coordinate(size_t j) const override;
  void print_info() const override;
};

class CircularDetector : public Detector_2D {
private:
  double distance;
  double R_min, R_max;
  double dR, d_phi;

public:
  // Planar annulus at z_local = D in the detector's own local frame: x_local = r*cos(phi), y_local = r*sin(phi),
  // r in [R_min, R_max] (N_R points), phi in [0, 2*pi) (N_phi points), then rotated like the other detectors.
  CircularDetector(size_t N_R, size_t N_phi, double D, double r_min, double r_max, double dir_x, double dir_y,
                   double dir_z, const MathUtils::RealFourTensor& laser_rotation);

  std::string get_type_name() const override { return "CircularDetector"; }
  double get_row_coordinate(size_t i) const override;
  double get_col_coordinate(size_t j) const override;
  void print_info() const override;
};

}  // namespace Core::Detector