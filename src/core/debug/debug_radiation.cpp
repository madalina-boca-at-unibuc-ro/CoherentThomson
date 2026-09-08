#include "../include/debug/debug_radiation.hpp"

#include <array>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace Core::Debug {

namespace {

// Mirrors Radiation::compute_radiation's internal bivector_element (radiation.cpp) -- deliberately
// duplicated rather than shared, so this debug path stays fully isolated from the production one.
inline double bivector_element(const MathUtils::RealFourVector& n, const MathUtils::RealFourVector& u, size_t alpha,
                               size_t beta) {
  return n[alpha] * u[beta] - n[beta] * u[alpha];
}

// The fixed (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) index order used throughout the codebase for the 6
// independent upper-triangle Faraday bivector elements (see MathUtils::ComplexBivector).
constexpr std::array<std::pair<size_t, size_t>, 6> kBivectorIndex = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};

}  // namespace

void export_radiation_integrand(const Particle::Electron& electron, const MathUtils::RealFourVector& detector_point,
                                double k, const std::string& filepath) {
  const std::vector<Particle::Electron::State>& trajectory = electron.get_trajectory();
  if (trajectory.empty()) {
    throw std::runtime_error("Debug::export_radiation_integrand: electron has no recorded trajectory -- call "
                             "electron.compute_trajectory() before exporting the debug integrand");
  }

  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Debug::export_radiation_integrand: failed to open file for export: " + filepath);
  }

  file << std::scientific << std::setprecision(6);
  file << "# debug radiation integrand: one row per trajectory point (tau), for a single electron/screen "
          "point/wavenumber k (= omega/c)\n";
  file << "# per-tau contribution to the packed Faraday bivector, (n_R0^alpha u^beta - n_R0^beta u^alpha) * "
          "amp_{long,short}(tau) -- the raw terms Radiation::compute_radiation sums over tau to build the final "
          "field\n";
  file << "# LR/SR = long_range/short_range, printed as 're im' pairs, in the fixed (01,02,03,12,13,23) bivector "
          "order\n";
  file << "i_tau tau";
  for (auto [alpha, beta] : kBivectorIndex)
    file << " LR_F" << alpha << beta << "_re LR_F" << alpha << beta << "_im";
  for (auto [alpha, beta] : kBivectorIndex)
    file << " SR_F" << alpha << beta << "_re SR_F" << alpha << beta << "_im";
  file << "\n";

  for (size_t i_tau = 0; i_tau < trajectory.size(); ++i_tau) {
    const MathUtils::RealFourVector& x = trajectory[i_tau].position;
    const MathUtils::RealFourVector& u = trajectory[i_tau].momentum;

    // Same n_R0/R/amp_long/amp_short construction as Radiation::compute_radiation's inner tau loop --
    // see radiation.cpp for the derivation/comments. n_R0 = R_0/|R_0| = (x_0 - r_0)/|x_0 - r_0|
    // points from the emitting particle to the observer (detector_point - x, not x - detector_point).
    MathUtils::RealFourVector n_R0 = detector_point - x;
    double R = MathUtils::create_unit_light_like_vector_in_place(n_R0);

    double n_R0_contract_u = MathUtils::contract(n_R0, u);
    double n_R0_dot3_u = MathUtils::dot3(n_R0, u);
    double amp_short_0 = n_R0_dot3_u / (R * R * n_R0_contract_u);

    double phase_base = x[0] + R;
    MathUtils::Complex cexp = std::polar(1.0, phase_base * k);
    MathUtils::Complex amp_long = MathUtils::Complex{0.0, -k / R} * cexp;
    MathUtils::Complex amp_short = amp_short_0 * cexp;

    file << i_tau << " " << trajectory[i_tau].tau;
    for (auto [alpha, beta] : kBivectorIndex) {
      MathUtils::Complex long_term = amp_long * bivector_element(n_R0, u, alpha, beta);
      file << " " << long_term.real() << " " << long_term.imag();
    }
    for (auto [alpha, beta] : kBivectorIndex) {
      MathUtils::Complex short_term = amp_short * bivector_element(n_R0, u, alpha, beta);
      file << " " << short_term.real() << " " << short_term.imag();
    }
    file << "\n";
  }

  std::cout << "Successfully exported debug radiation integrand to " << filepath << "\n";
}

void export_radiation_phase(const Particle::Electron& electron, const MathUtils::RealFourVector& detector_point,
                            double k, const std::string& filepath) {
  const std::vector<Particle::Electron::State>& trajectory = electron.get_trajectory();
  if (trajectory.empty()) {
    throw std::runtime_error("Debug::export_radiation_phase: electron has no recorded trajectory -- call "
                             "electron.compute_trajectory() before exporting the debug phase");
  }

  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Debug::export_radiation_phase: failed to open file for export: " + filepath);
  }

  file << std::scientific << std::setprecision(6);
  file << "# debug radiation phase: one row per trajectory point (tau), for a single electron/screen "
          "point/wavenumber k (= omega/c)\n";
  file << "# phase = (x[0] + R) * k, unwrapped (not reduced mod 2*pi); exp_re/exp_im = cos(phase)/sin(phase) -- "
          "the exp(i*phase) factor common to every long-range/short-range bivector component in "
          "debug_integrand.dat\n";
  file << "i_tau tau phase exp_re exp_im\n";

  for (size_t i_tau = 0; i_tau < trajectory.size(); ++i_tau) {
    const MathUtils::RealFourVector& x = trajectory[i_tau].position;

    MathUtils::RealFourVector n_R0 = detector_point - x;
    double R = MathUtils::create_unit_light_like_vector_in_place(n_R0);

    double phase_base = x[0] + R;
    double phase = phase_base * k;
    MathUtils::Complex cexp = std::polar(1.0, phase);

    file << i_tau << " " << trajectory[i_tau].tau << " " << phase << " " << cexp.real() << " " << cexp.imag() << "\n";
  }

  std::cout << "Successfully exported debug radiation phase to " << filepath << "\n";
}

}  // namespace Core::Debug
