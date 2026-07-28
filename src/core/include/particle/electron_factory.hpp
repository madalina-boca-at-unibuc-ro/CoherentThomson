#pragma once
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../io_utils/io_utils.hpp"
#include "../math_utils/math_constants.hpp"
#include "../phys_utils/phys_utils.hpp"
#include "electron.hpp"

namespace Core::Particle {

// Cylinder-beam config parsing lives in IoUtils (shared with
// Simulation::init_simulation_parameters, which needs the same average_px/py/pz).
using CylinderBeamParams = Core::IoUtils::CylinderBeamParams;

// Draws one electron: position uniform over the cylinder cross-section
// (area-preserving in r) and uniform along the cylinder's height, momentum
// Gaussian about params.average_p*. The cylinder's height axis is first
// reoriented via beam_axis_rotation to match the beam's own mean momentum
// direction (so the finite-height extent lies along the direction the beam
// actually moves in, rather than always along canonical Oz regardless of
// where average_p* points), then both position and momentum are rotated
// into the lab frame via rotation_matrix (the laser's own rotation, shared
// with the beam). Advances gen's state, so a caller generating a whole beam
// must share one generator across calls to get a single reproducible
// sequence.
Electron generate_electron(const CylinderBeamParams& params, const Core::MathUtils::RealFourTensor& rotation_matrix,
                           const Core::MathUtils::RotationMatrix3x3& beam_axis_rotation, double tau_0, double d_tau,
                           size_t N_tau, std::mt19937& gen, bool store_trajectory);

std::vector<Electron> generate_cylinder_beam(const ConfigMap& config,
                                             const Core::MathUtils::RealFourTensor rotation_matrix, const double tau_0,
                                             const double d_tau, const size_t N_tau);

}  // namespace Core::Particle
