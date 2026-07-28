#pragma once
#include <tuple>

#include "../math_utils/math_utils.hpp"

namespace Core::Laser {

// Clean alias for a 4x4 matrix tensor representation
using FaradayTensor = Core::MathUtils::FourTensor<double>;

class LaserField {
protected:
  // Laser parameters
  double omega;          // Laser angular frequency
  double a0;             // Dimensionless laser strength parameter (a0 = eE0/(m_e c omega))
  double E0_c;           // Peak electric field amplitude divided by c for proper unit
                         // consistency in the Faraday tensor
  double flat_duration;  // Half-width of the flat-top plateau phase
  double wing_sigma, delay,
      wing_sigma_cutoff;              // Standard deviation phase width of Gaussian wings
                                      // the overall delay shift in phase units and the cutoff factor of the Gaussian
  MathUtils::Complex zeta_1, zeta_2;  // Complex polarization coefficients for the two orthogonal
                                      // polarization directions; e.g. zeta_1=(1,0), zeta_2=(0,1) gives
                                      // circular polarization

  // Parameters to be used in the simulation
  double phi_min, phi_max;  // Phase range for the laser pulse and any global
                            // delay shift in phase units
  double d_phi;             // phase step used in laser plotting only
  size_t NT, N_phi;         // number of points per cycle / total number of points used in laser plotting only
                            // NT is read from the .cfg file d_phi = 2 * pi /NT, N_phi = (phi_max-phi_min)/d_phi

  Core::MathUtils::RealFourVector unity_n;               // Normalized propagation direction 4-vector (n)
  Core::MathUtils::RealFourVector epsilon_1, epsilon_2;  // Normalized polarization direction vectors (\hat{\epsilon});
                                                         // they must be orthogonal on k
  Core::MathUtils::RealFourTensor rotation_matrix;

  // Internal helper to calculate the temporal envelope factor g(\phi); shared by every derived laser
  // type -- the Gaussian-flat-top pulse shape in time is the same regardless of the transverse mode.
  // returns only the exponant of the envelope, will be combined in the calcualtion of the amplitude
  double envelope(double phi) const;

public:
  // Constructor initializing the baseline laser characteristics shared by every wave type: temporal
  // envelope, propagation direction, polarization, and the resulting rotation matrices.
  LaserField(double angular_freq, double a0, double pulse_flat_phase, double wing_sigma_phase, double delay_in,
             double wing_cutoff_sigmas, MathUtils::Complex zeta_1_in, MathUtils::Complex zeta_2_in, size_t NT_in,
             const Core::MathUtils::RealFourVector& n_dir_in);
  virtual ~LaserField() = default;

  // Computes and returns the complete F^{\mu\nu} tensor at a given 4-vector position. This is the one
  // piece of physics that differs between wave types (plane wave vs. Laguerre-Gauss, ...).
  virtual FaradayTensor get_faraday_tensor(const Core::MathUtils::RealFourVector& x_mu) const = 0;

  double get_omega() const { return omega; }
  double get_a0() const { return a0; }
  double get_flat_duration() const { return flat_duration; }
  double get_wing_sigma_cutoff() const { return wing_sigma_cutoff; }
  MathUtils::Complex get_zeta_1() const { return zeta_1; }
  MathUtils::Complex get_zeta_2() const { return zeta_2; }
  double get_phi_min() const { return phi_min; }
  double get_phi_max() const { return phi_max; }
  double get_d_phi() const { return d_phi; }
  size_t get_NT() const { return NT; }
  size_t get_N_phi() const { return N_phi; }
  const Core::MathUtils::RealFourVector get_unity_n() const { return unity_n; }
  const Core::MathUtils::RealFourVector get_epsilon_1() const { return epsilon_1; }
  const Core::MathUtils::RealFourVector get_epsilon_2() const { return epsilon_2; }
  const Core::MathUtils::RealFourTensor get_rotation_matrix() const { return rotation_matrix; }
};

// The original plane-wave laser: a spatially uniform transverse profile modulated only by the temporal
// envelope. Reuses the base constructor as-is (no additional parameters).
class PlaneWaveLaser : public LaserField {
private:
  MathUtils::Complex complex_amplitude(const MathUtils::RealFourVector& x) const;

public:
  using LaserField::LaserField;
  FaradayTensor get_faraday_tensor(const Core::MathUtils::RealFourVector& x_mu) const override;
};

// A Laguerre-Gauss (LG_{p,l}) laser mode: adds a transverse radial/azimuthal profile (with Gouy phase
// and wavefront curvature) on top of the same temporal envelope and carrier phase used by
// PlaneWaveLaser.
class LaguerreGaussLaser : public LaserField {
private:
  int p;       // Radial index (number of radial nodes)
  int l;       // Azimuthal (orbital angular momentum) index
  double w0;   // Beam waist at the focus
  double z_R;  // Rayleigh range, derived from w0 and omega: z_R = pi * w0^2 / lambda
  // Returns {amplitude, d(amplitude)/dx_loc, d(amplitude)/dy_loc} at x_mu, all evaluated together since
  // they share almost all of their computation; x_loc/y_loc are the laser's canonical-frame transverse
  // coordinates. The x/y derivatives are needed (unlike in the plane-wave case) to build the E and B
  // field components along the propagation direction from the divergence-free condition.
  std::tuple<MathUtils::Complex, MathUtils::Complex, MathUtils::Complex> complex_amplitude(
      const MathUtils::RealFourVector& x) const;

public:
  LaguerreGaussLaser(double angular_freq, double a0, double pulse_flat_phase, double wing_sigma_phase, double delay_in,
                     double wing_cutoff_sigmas, MathUtils::Complex zeta_1_in, MathUtils::Complex zeta_2_in,
                     size_t NT_in, const Core::MathUtils::RealFourVector& n_dir_in, int p_in, int l_in, double w0_in);
  FaradayTensor get_faraday_tensor(const Core::MathUtils::RealFourVector& x_mu) const override;

  int get_lg_p() const { return p; }
  int get_lg_l() const { return l; }
  double get_lg_w0() const { return w0; }
  double get_lg_rayleigh_range() const { return z_R; }
};

}  // namespace Core::Laser
