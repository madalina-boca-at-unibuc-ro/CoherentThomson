#pragma once

#include <string>
#include <vector>

#include "../simulation/simulation.hpp"

namespace Core::Radiation {

// Exports every element of the coherently-summed long_range/short_range Faraday tensors, one row per
// (frequency, screen-point) pair, for later inspection/plotting.
void plot_radiation_field(const Simulation::RadiationField& field, const std::vector<double>& frequencies_list,
                          const std::string& filepath);

}  // namespace Core::Radiation
