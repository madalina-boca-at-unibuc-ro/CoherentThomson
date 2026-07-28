#pragma once

#include "../math_utils/math_utils.hpp"
#include "../particle/electron.hpp"
#include "../simulation/simulation.hpp"

namespace Core::Radiation {

// Computes one electron's contribution to the radiated field and adds it into `field` (on top of
// whatever it already holds) at every configured frequency/screen-point pair. `field` holds the
// packed accumulation-time representation (Simulation::PackedFaraday); run_simulation
// reconstructs the true Faraday tensor (Simulation::Faraday) once every electron's and every
// thread's contribution has been summed.
void compute_radiation(Particle::Electron& electron, const Laser::LaserField& laser,
                       const std::vector<double>& frequencies_list, const Detector::Detector_2D& detector,
                       Simulation::PackedRadiationField& field);

}  // namespace Core::Radiation
