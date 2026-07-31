#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../core/include/detector/detector_factory.hpp"
#include "../core/include/detector/detector_plotter.hpp"
#include "../core/include/io_utils/config_parser.hpp"
#include "../core/include/io_utils/io_utils.hpp"
#include "../core/include/laser/laser_factory.hpp"
#include "../core/include/laser/laser_plotter.hpp"
#include "../core/include/particle/electron_factory.hpp"
#include "../core/include/particle/electron_plotter.hpp"
#include "../core/include/radiation/radiation_plotter.hpp"
#include "../core/include/simulation/simulation.hpp"

int main(int argc, char* argv[]) {
  using namespace Core;
  if (argc < 2) {
    std::cerr << "Error: No configuration file provided.\n";
    std::cerr << "Usage: " << argv[0] << " <config_file>\n";
    return 1;
  }

  std::string config_path = argv[1];
  try {
    // Read the configuration file and print the configurations parsed.
    ConfigMap simulation_config = IoUtils::read_config_file(config_path, true);

    std::string run_output_dir = IoUtils::make_run_output_directory(simulation_config);
    std::cout << "Writing run output to " << run_output_dir << "\n";
    IoUtils::copy_config_to_run_directory(config_path, run_output_dir);

    auto laser = Laser::create_laser(simulation_config);
    Laser::export_field_vs_phase(*laser, run_output_dir + "/laser_field.dat");

    if (IoUtils::get_required(simulation_config, "plot_field_heatmap") == "true") {
      std::string heatmap_axes_unit = "lambda";
      auto heatmap_axes_scale = IoUtils::convert_unit_to_number(heatmap_axes_unit, simulation_config);

      // A plane wave has no transverse profile at all (the heat map is uniform regardless of window
      // size), so its window stays fixed at the field_heatmap_x/y_min/max config values below. A
      // Laguerre-Gauss mode's actual transverse scale is set by its own waist w0, not by those fixed
      // values, so size the window to it instead (+-2*w0, generous enough to show the mode's rings/
      // lobes for the small p/l indices this is meant to be used with).
      std::string laser_type = IoUtils::get_required(simulation_config, "laser_type");
      double x_min, x_max, y_min, y_max;
      if (laser_type == "laguerre_gauss") {
        double lg_w0 = std::get<2>(IoUtils::get_laser_lg_params(simulation_config));
        x_min = -2.0 * lg_w0;
        x_max = 2.0 * lg_w0;
        y_min = -2.0 * lg_w0;
        y_max = 2.0 * lg_w0;
      } else {
        auto [x_min_val, x_min_unit] =
            IoUtils::split_value_and_unit(IoUtils::get_required(simulation_config, "field_heatmap_x_min"));
        x_min = x_min_val * IoUtils::convert_unit_to_number(x_min_unit, simulation_config);
        auto [x_max_val, x_max_unit] =
            IoUtils::split_value_and_unit(IoUtils::get_required(simulation_config, "field_heatmap_x_max"));
        x_max = x_max_val * IoUtils::convert_unit_to_number(x_max_unit, simulation_config);
        auto [y_min_val, y_min_unit] =
            IoUtils::split_value_and_unit(IoUtils::get_required(simulation_config, "field_heatmap_y_min"));
        y_min = y_min_val * IoUtils::convert_unit_to_number(y_min_unit, simulation_config);
        auto [y_max_val, y_max_unit] =
            IoUtils::split_value_and_unit(IoUtils::get_required(simulation_config, "field_heatmap_y_max"));
        y_max = y_max_val * IoUtils::convert_unit_to_number(y_max_unit, simulation_config);
      }
      size_t heatmap_Nx = std::stoull(IoUtils::get_required(simulation_config, "field_heatmap_Nx"));
      size_t heatmap_Ny = std::stoull(IoUtils::get_required(simulation_config, "field_heatmap_Ny"));

      // The actual instant the flat-top plateau begins at: LaserField's constructor (laser_field.cpp)
      // sets its internal phase-domain `delay` to wing_sigma_cutoff * wing_sigma (currently ignoring
      // the parsed laser_delay entirely -- see CLAUDE.md), and envelope(phi) switches from the leading
      // Gaussian wing to the flat plateau exactly at phi = delay. Recomputing that same product here
      // (rather than using laser_delay) keeps this in sync with the real switch point regardless of
      // whether laser_delay happens to match wing_sigma_cutoff * wing_sigma.
      double t_heatmap = IoUtils::get_laser_wing_sigma_cutoff(simulation_config) *
                         IoUtils::get_laser_wing_sigma(simulation_config) / laser->get_omega();

      Laser::export_field_heatmap_z0(*laser, t_heatmap, x_min, x_max, y_min, y_max, heatmap_Nx, heatmap_Ny,
                                     heatmap_axes_scale, heatmap_axes_unit,
                                     run_output_dir + "/laser_field_heatmap_z0.dat");
      std::cout << "Successfully exported laser field heatmap to " << run_output_dir << "/laser_field_heatmap_z0.dat\n";
    }

    Simulation::simulation_parameters sim_par = Simulation::init_simulation_parameters(simulation_config, *laser);

    auto electron_beam =
        Particle::generate_cylinder_beam(simulation_config, laser->get_rotation_matrix(), sim_par.tau_0_traj,
                                         sim_par.d_tau_traj, sim_par.simulation_length);

    if (IoUtils::get_required(simulation_config, "plot_beam_scatter") == "true") {
      std::string beam_axes_unit = "lambda";
      auto beam_axes_scale = IoUtils::convert_unit_to_number(beam_axes_unit, simulation_config);
      Particle::plot_electron_beam_scatter(electron_beam, run_output_dir + "/electron_beam_scatter.dat",
                                           beam_axes_scale, beam_axes_unit);
    }

    size_t num_trajectory_electrons = std::min<size_t>(10, electron_beam.size());
    std::vector<size_t> trajectory_electron_indices(num_trajectory_electrons);
    for (size_t i = 0; i < num_trajectory_electrons; ++i) {
      trajectory_electron_indices[i] = i;
      electron_beam[i].compute_trajectory(*laser);
    }

    Particle::plot_particle_trajectory(electron_beam, trajectory_electron_indices, run_output_dir + "/electron.dat");

    auto detector = Detector::create_detector(simulation_config, *laser);
    std::string detector_axes_unit = "lambda";
    auto detector_axes_scale = IoUtils::convert_unit_to_number(detector_axes_unit, simulation_config);
    Detector::plot_detector(*detector, "stereographic_plot", run_output_dir, detector_axes_scale, detector_axes_unit);
    if (IoUtils::get_required(simulation_config, "plot_detector_scatter") == "true") {
      Detector::plot_detector_scatter(*detector, run_output_dir, detector_axes_scale, detector_axes_unit);
    }

    if (IoUtils::get_required(simulation_config, "dense_frequency_spectrum") == "true" &&
        detector->get_total_points() != 1) {
      std::cerr << "Warning: dense_frequency_spectrum=true is meant for a detector collapsed to a single "
                << "point, but this detector has " << detector->get_total_points() << " points -- the fine "
                << "frequency scan will run over every one of them.\n";
    }

    size_t num_threads = IoUtils::get_num_threads(simulation_config);
    auto simulation_start = std::chrono::steady_clock::now();
    Simulation::RadiationField radiation_field = Simulation::run_simulation(
        simulation_config, *laser, *detector, electron_beam, sim_par.frequencies, num_threads);
    auto simulation_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> simulation_elapsed = simulation_end - simulation_start;
    std::cout << "Simulation time:      " << simulation_elapsed.count() << " s\n";
    std::cout << "Number of particles:  " << electron_beam.size() << "\n";
    std::cout << "Number of screen pts: " << detector->get_total_points() << "\n";
    std::cout << "Number of threads:    " << num_threads << "\n";

    Radiation::plot_radiation_field(radiation_field, sim_par.frequencies, run_output_dir + "/radiation_field.dat");

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}