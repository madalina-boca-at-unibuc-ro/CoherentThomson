#pragma once
#include <cstddef>
#include <string>

#include "../detector/detector.hpp"
#include "../io_utils/config_map.hpp"
#include "../laser/laser_field.hpp"
#include "../simulation/simulation.hpp"

namespace Core::Logging {

// Writes a human-readable summary of this run's full parameter set -- both the raw config values and
// everything derived/computed from them (dressed momentum q, the actual numeric frequencies_list,
// detector geometry re-expressed in lambda/w0 units, ...) -- to `filepath`. Meant as a richer,
// at-a-glance companion to the verbatim config.cfg snapshot IoUtils::copy_config_to_run_directory
// already writes into every run directory: config.cfg has everything needed to *reproduce* a run,
// this has everything needed to *understand what was actually computed* without re-deriving it by
// hand. Re-reads geometry/beam keys straight from `config` (the same way detector_factory.cpp/
// electron_factory.cpp already do) rather than adding new getters to Detector/Electron, so this stays
// purely additive and does not touch the physics classes themselves. `num_threads` should be the
// value Simulation::run_simulation resolved (not the raw, possibly-zero-meaning-"all" config value),
// `simulation_elapsed_seconds` the wall-clock time that run took, and `total_cpu_seconds` the same
// summed-across-threads CPU-seconds figure Simulation::run_simulation resolved (its own out-param) --
// used here (not simulation_elapsed_seconds) to derive the per-electron/per-screen-point/per-cycle
// timing figure, for the same reason run_simulation itself prefers it over wall-clock time: dividing
// wall-clock time by electron/screen-point count would understate the true per-unit cost by roughly
// num_threads, since electrons are processed in parallel, not sequentially.
void write_run_log(const ConfigMap& config, const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                   size_t num_electrons, const Simulation::simulation_parameters& sim_par, size_t num_threads,
                   double simulation_elapsed_seconds, double total_cpu_seconds, const std::string& filepath);

}  // namespace Core::Logging
