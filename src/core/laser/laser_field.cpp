#include "../include/laser/laser_field.hpp"

#include "../include/math_utils/math_utils.hpp"
#include "../include/phys_utils/phys_utils.hpp"

namespace Core::Laser {

// =========================================================================
// BASE CLASS LASERFIELD DEFINITIONS
// =========================================================================

LaserField::LaserField(double omega_in, double a0_in, double flat_duration_in, double wing_sigma_in, double delay_in,
                       double wing_sigma_cutoff_in, MathUtils::Complex zeta_1_in, MathUtils::Complex zeta_2_in,
                       size_t NT_in, const Core::MathUtils::RealFourVector& unity_n_in)
    : omega(omega_in),
      a0(a0_in),
      flat_duration(flat_duration_in),
      wing_sigma(wing_sigma_in),
      delay(delay_in),
      wing_sigma_cutoff(wing_sigma_cutoff_in),
      zeta_1(zeta_1_in),
      zeta_2(zeta_2_in),
      NT(NT_in) {
  E0_c = a0 * omega * Core::PhysUtils::AtomicUnits::m_0 /
         (Core::PhysUtils::AtomicUnits::e_0);  // E0 = omega * a0 * m_e  / |e|
  // Create the normalized 4-vector for the propagation direction (with time
  // component fixed at 1.0 for normalization)
  unity_n = Core::MathUtils::create_unit_light_like_vector(unity_n_in);

  // Create the polarization vectors (orthogonal to n) and normalize them. For
  // simplicity we start with two polarization vecotrs along Ox, Oy, and the
  // apply the rotation with the polar angles of unity_n to ensure they are
  // orthogonal to the propagation direction.

  Core::MathUtils::RealFourVector pol_dir1(0.0, 1.0, 0.0, 0.0);
  Core::MathUtils::RealFourVector pol_dir2(0.0, 0.0, 1.0, 0.0);

  // Apply rotation to ensure polarization vectors are orthogonal to the
  // propagation direction
  rotation_matrix = Core::MathUtils::rotation_four_tensor_from_direction(unity_n[1], unity_n[2], unity_n[3]);
  epsilon_1 = Core::MathUtils::contract(rotation_matrix, pol_dir1);
  epsilon_2 = Core::MathUtils::contract(rotation_matrix, pol_dir2);
  delay = wing_sigma_cutoff *
          wing_sigma;  // set delay such that the simulation starts form zero.
                       // TO DO: To improve the choice, to make is consisten with the initial simulation time
  phi_min = -wing_sigma_cutoff * wing_sigma + delay;                 // Start well in advance of the leading wing
  phi_max = flat_duration + wing_sigma_cutoff * wing_sigma + delay;  // End well after the trailing wing
  d_phi = 2 * MathUtils::pi / NT;
  N_phi = static_cast<size_t>(
      std::round((phi_max - phi_min) / d_phi));  // Number of phi points to cover the full range with NT points per
}

// returns only the explonent of the envelope, to be used in the calculation of the complex amplitude
double LaserField::envelope(double phi) const {
  if (phi < delay) {
    // Leading Gaussian wing
    double delta_phi = phi - delay;
    return (-(delta_phi * delta_phi) / (2.0 * wing_sigma * wing_sigma));
  } else if (phi > flat_duration + delay) {
    // Trailing Gaussian wing
    double delta_phi = phi - (flat_duration + delay);
    return (-(delta_phi * delta_phi) / (2.0 * wing_sigma * wing_sigma));
  } else {
    // Constant flat-top plateau
    return 0.0;
  }
}

// =========================================================================
// PLANE WAVE LASER IMPLEMENTATION
// =========================================================================

MathUtils::Complex PlaneWaveLaser::complex_amplitude(const MathUtils::RealFourVector& x) const {
  double phi = omega / PhysUtils::AtomicUnits::c *
               MathUtils::contract(x, unity_n);  // Using the dot product of x_mu and n to get the phase
  MathUtils::Complex exponent = envelope(phi) + MathUtils::I * phi;
  MathUtils::Complex result = E0_c * std::exp(exponent);
  return result;
}

FaradayTensor PlaneWaveLaser::get_faraday_tensor(const Core::MathUtils::RealFourVector& x_mu) const {
  MathUtils::Complex amplitude = complex_amplitude(x_mu);
  // E = Real(epsilon_1 * zeta_1 * amplitude + epsilon_2 * zeta_2 * amplitude); zeta_1/zeta_2 complex so
  // e.g. zeta_1=(1,0), zeta_2=(0,1) gives circular polarization.
  MathUtils::Complex weighted_amplitude_1 = zeta_1 * amplitude;
  MathUtils::Complex weighted_amplitude_2 = zeta_2 * amplitude;

  // Electric Field Vector Components
  double Ex_c = epsilon_1[1] * real(weighted_amplitude_1) + epsilon_2[1] * real(weighted_amplitude_2);
  double Ey_c = epsilon_1[2] * real(weighted_amplitude_1) + epsilon_2[2] * real(weighted_amplitude_2);
  double Ez_c = epsilon_1[3] * real(weighted_amplitude_1) + epsilon_2[3] * real(weighted_amplitude_2);

  // Magnetic Field Vector Components via cross product: B = (k_dir \times E) /
  // c
  double Bx = (unity_n[2] * Ez_c - unity_n[3] * Ey_c);
  double By = (unity_n[3] * Ex_c - unity_n[1] * Ez_c);
  double Bz = (unity_n[1] * Ey_c - unity_n[2] * Ex_c);

  // 4. Construct the antisymmetric Faraday Tensor matrix F^{\mu\nu}
  FaradayTensor F;

  // Row 0 (t)
  F[0][0] = 0.0;
  F[0][1] = -Ex_c;
  F[0][2] = -Ey_c;
  F[0][3] = -Ez_c;
  // Row 1 (x)
  F[1][0] = Ex_c;
  F[1][1] = 0.0;
  F[1][2] = -Bz;
  F[1][3] = By;
  // Row 2 (y)
  F[2][0] = Ey_c;
  F[2][1] = Bz;
  F[2][2] = 0.0;
  F[2][3] = -Bx;
  // Row 3 (z)
  F[3][0] = Ez_c;
  F[3][1] = -By;
  F[3][2] = Bx;
  F[3][3] = 0.0;

  return F;
}

// =========================================================================
// LAGUERRE-GAUSS LASER IMPLEMENTATION
// =========================================================================

LaguerreGaussLaser::LaguerreGaussLaser(double angular_freq, double a0_in, double pulse_flat_phase,
                                       double wing_sigma_phase, double delay_in, double wing_cutoff_sigmas,
                                       MathUtils::Complex zeta_1_in, MathUtils::Complex zeta_2_in, size_t NT_in,
                                       const Core::MathUtils::RealFourVector& n_dir_in, int p_in, int l_in,
                                       double w0_in)
    : LaserField(angular_freq, a0_in, pulse_flat_phase, wing_sigma_phase, delay_in, wing_cutoff_sigmas, zeta_1_in,
                 zeta_2_in, NT_in, n_dir_in),
      p(p_in),
      l(l_in),
      w0(w0_in) {
  // Rayleigh range z_R = pi * w0^2 / lambda = w0^2 * omega / (2c), consistent with lambda = 2*pi*c/omega.
  z_R = (w0 * w0 * omega) / (2.0 * Core::PhysUtils::AtomicUnits::c);
}

// Returns {amplitude, d(amplitude)/dx_loc, d(amplitude)/dy_loc} together, since they share almost all
// of their computation. Unlike the plane-wave case, the field components along the propagation
// direction need these transverse derivatives to enforce div(E) = 0, div(B) = 0.
//
// The azimuthal factor rho^|l| * exp(i*l*azimuth) is built as the exact complex polynomial
// (x_loc + i*sign(l)*y_loc)^|l|, identically equal to it since rho*exp(i*sign(l)*azimuth) = x_loc +
// i*sign(l)*y_loc -- this stays smooth at rho=0 and lets the x/y derivatives be plain polynomial
// derivatives rather than polar-coordinate expressions that are singular on-axis.
std::tuple<MathUtils::Complex, MathUtils::Complex, MathUtils::Complex> LaguerreGaussLaser::complex_amplitude(
    const MathUtils::RealFourVector& x_mu) const {
  // Carrier + temporal-envelope phase -- identical construction to the plane-wave case; the Minkowski
  // contraction is invariant under rotating x_mu and unity_n together, so this is exactly the on-axis
  // phase phi = omega*t - k*z regardless of frame.
  double phi = omega / Core::PhysUtils::AtomicUnits::c * Core::MathUtils::contract(x_mu, unity_n);

  // Transverse position and longitudinal distance from the beam waist, in the laser's own canonical
  // frame (propagation along local +Z); the waist is assumed to sit at the canonical-frame origin.
  // epsilon_1/epsilon_2/unity_n are the canonical frame's x/y/z axes expressed in the lab frame (the
  // same columns that used to make up the now-removed spatial_rotation matrix), so dotting x_mu's
  // spatial part against each is equivalent to rotating x_mu into that frame.
  double x_loc = Core::MathUtils::dot3(x_mu, epsilon_1);
  double y_loc = Core::MathUtils::dot3(x_mu, epsilon_2);
  double z_loc = Core::MathUtils::dot3(x_mu, unity_n);
  double s = x_loc * x_loc + y_loc * y_loc;  // = rho^2
  double k_wave = omega / Core::PhysUtils::AtomicUnits::c;

  int n = std::abs(l);
  double Npn = std::sqrt(2.0) / MathUtils::factorial(n) *
               std::sqrt(MathUtils::factorial(p + n) / MathUtils::factorial(p));  // the normalization constant
  double w_z = w0 * std::sqrt(1.0 + (z_loc / z_R) * (z_loc / z_R));               // the beam radius at z
  double C_n = E0_c * Npn * (w0 / w_z) * std::pow(std::sqrt(2.0) / w_z, n);       // prefactor
  double u = 2.0 * s / (w_z * w_z);  // argument of the hypergeometric function
  double hypergeometric_val = MathUtils::hypergeometric_1F1_neg_int_a(-p, n + 1, u);  // the hypergeomatric function

  // Vortex factor V = (x_loc + i*sign(l)*y_loc)^|l|, identically equal to rho^|l| * exp(i*l*azimuth).
  double sign_l = (l >= 0) ? 1.0 : -1.0;
  MathUtils::Complex zeta_c(x_loc, sign_l * y_loc);
  MathUtils::Complex V = std::pow(zeta_c, n);

  // the imaginary part of the phase
  double gouy_phase = -static_cast<double>(2 * p + n + 1) * std::atan(z_loc / z_R);
  double kappa_phase = k_wave * s * z_loc /
                       (2.0 * (z_loc * z_loc + z_R * z_R));  // curvature_phase = kappa * z * rho^2 / (2 ( z^2 + z_R^2))

  double gaussian_phase = -s / (w_z * w_z);
  double envelope_phase = envelope(phi);
  MathUtils::Complex prefactor =
      std::exp(MathUtils::Complex(envelope_phase + gaussian_phase, phi + gouy_phase + kappa_phase));

  MathUtils::Complex amplitude = prefactor * V * C_n * hypergeometric_val;
  MathUtils::Complex d_amplitude_dx = 0.0;
  MathUtils::Complex d_amplitude_dy = 0.0;

  return {amplitude, d_amplitude_dx, d_amplitude_dy};
}

// The main difference with respect to the plane wave case is the presence of field components along the propagation
// direction. They are chosen such that div(E) = 0; div(B) = 0.

// The complex scalar solution is built in the complex_amplitude

FaradayTensor LaguerreGaussLaser::get_faraday_tensor(const Core::MathUtils::RealFourVector& x_mu) const {
  MathUtils::Complex amplitude = std::get<0>(complex_amplitude(x_mu));
  // Same polarization construction as PlaneWaveLaser (see there); this is currently a stub identical
  // to the plane wave, pending the true LG transverse-mode physics.
  MathUtils::Complex weighted_amplitude_1 = zeta_1 * amplitude;
  MathUtils::Complex weighted_amplitude_2 = zeta_2 * amplitude;

  // Electric Field Vector Components
  double Ex_c = epsilon_1[1] * real(weighted_amplitude_1) + epsilon_2[1] * real(weighted_amplitude_2);
  double Ey_c = epsilon_1[2] * real(weighted_amplitude_1) + epsilon_2[2] * real(weighted_amplitude_2);
  double Ez_c = epsilon_1[3] * real(weighted_amplitude_1) + epsilon_2[3] * real(weighted_amplitude_2);

  // Magnetic Field Vector Components via cross product: B = (k_dir \times E) /
  // c
  double Bx = (unity_n[2] * Ez_c - unity_n[3] * Ey_c);
  double By = (unity_n[3] * Ex_c - unity_n[1] * Ez_c);
  double Bz = (unity_n[1] * Ey_c - unity_n[2] * Ex_c);

  // 4. Construct the antisymmetric Faraday Tensor matrix F^{\mu\nu}
  FaradayTensor F;

  // Row 0 (t)
  F[0][0] = 0.0;
  F[0][1] = -Ex_c;
  F[0][2] = -Ey_c;
  F[0][3] = -Ez_c;
  // Row 1 (x)
  F[1][0] = Ex_c;
  F[1][1] = 0.0;
  F[1][2] = -Bz;
  F[1][3] = By;
  // Row 2 (y)
  F[2][0] = Ey_c;
  F[2][1] = Bz;
  F[2][2] = 0.0;
  F[2][3] = -Bx;
  // Row 3 (z)
  F[3][0] = Ez_c;
  F[3][1] = -By;
  F[3][2] = Bx;
  F[3][3] = 0.0;

  return F;
}

}  // namespace Core::Laser
