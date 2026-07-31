#pragma once

#include <string>
#include <vector>

#include "../simulation/simulation.hpp"

namespace Core::Radiation {

// Exports every element of the coherently-summed long_range/short_range Faraday tensors, one row per
// (frequency, screen-point) pair, for later inspection/plotting. The exported "omega" column is
// frequencies_list normalized by fundamental_frequency (Simulation::simulation_parameters::fundamental_frequency,
// i.e. PhysUtils::non_linear_Thomson_formula(k1, p, n2, 1)) rather than left in raw atomic-unit omega, so a
// value of 3.0 reads as "third harmonic" regardless of the observation direction's Doppler shift.
void plot_radiation_field(const Simulation::RadiationField& field, const std::vector<double>& frequencies_list,
                          double fundamental_frequency, const std::string& filepath);

}  // namespace Core::Radiation
