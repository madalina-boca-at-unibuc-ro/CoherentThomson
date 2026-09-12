#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../core/include/debug/debug_radiation.hpp"
#include "../core/include/detector/detector_factory.hpp"
#include "../core/include/detector/detector_plotter.hpp"
#include "../core/include/io_utils/config_parser.hpp"
#include "../core/include/io_utils/io_utils.hpp"
#include "../core/include/laser/laser_factory.hpp"
#include "../core/include/laser/laser_plotter.hpp"
#include "../core/include/logging/run_log.hpp"
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

      // Size the window to the electron beam's own transverse extent (+-beam_cylinder_radius)
      // instead of a laser-type-specific heuristic -- the point of this heatmap is to see the field
      // where the beam actually sits, and the beam's radius is a laser-type-independent quantity
      // already available straight from the config (no need to wait for generate_cylinder_beam,
      // which hasn't run yet at this point in main -- IoUtils::parse_cylinder_beam_params reads the
      // same beam_cylinder_radius key it would use). Deliberately ignores beam_center_x/y: those are
      // applied in the beam's own local frame, before the mean-momentum rotation
      // generate_cylinder_beam performs, which generally does *not* land back on this heatmap's
      // canonical (laser-along-Oz) frame unless the beam's mean momentum has zero transverse
      // component -- centering the window would need reproducing that rotation here, so the window
      // stays beam-radius-sized but always centered on the axis.
      // Falls back to the previous windowing (fixed field_heatmap_x/y_min/max for plane_wave,
      // +-2*w0 for laguerre_gauss, whose actual transverse scale is set by its own waist rather than
      // any beam property) only when the beam has no meaningful transverse extent to size off of
      // (radius == 0, e.g. the single on-axis point beam in config/coherent_thomson_debug.cfg) --
      // a zero-radius window would otherwise collapse the heatmap to a single point.
      IoUtils::CylinderBeamParams beam_params = IoUtils::parse_cylinder_beam_params(simulation_config);
      double x_min, x_max, y_min, y_max;
      if (beam_params.radius > 0.0) {
        x_min = -beam_params.radius;
        x_max = beam_params.radius;
        y_min = -beam_params.radius;
        y_max = beam_params.radius;
      } else {
        std::string laser_type = IoUtils::get_required(simulation_config, "laser_type");
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

    // Debug mode: dump the per-tau long-range/short-range radiation integrand for a single
    // electron/screen point/frequency (the fundamental, i.e. N_harmonics=1 regardless of
    // dense_frequency_spectrum/N_harmonics/N_omega) -- see Core::Debug::export_radiation_integrand.
    // Purely additive: it only reads already-built state and does not alter the rest of the run, so
    // the normal pipeline below (including run_simulation) still executes and radiation_field.dat
    // can be cross-checked against the integrand's own tau-sum.
    if (IoUtils::get_required(simulation_config, "debug") == "true") {
      if (electron_beam.size() != 1) {
        throw std::runtime_error("debug=true requires exactly one electron in the beam (set beam_particle_count=1)");
      }
      if (detector->get_type_name() != "RectangularDetector") {
        throw std::runtime_error("debug=true requires a rectangular detector (set detector_type=rectangular)");
      }
      if (detector->get_total_points() != 1) {
        throw std::runtime_error("debug=true requires exactly one detector point (set rectangular_detector_Nx="
                                 "rectangular_detector_Ny=1)");
      }
      if (IoUtils::get_required(simulation_config, "radiation_formula") != "simplified") {
        throw std::runtime_error("debug=true requires radiation_formula=simplified: "
                                 "Debug::export_radiation_integrand/export_radiation_phase only reimplement the "
                                 "simplified form's per-tau math, so they cannot cross-check a direct-formula run");
      }
      Debug::export_radiation_integrand(electron_beam[0], detector->get_point(0), sim_par.fundamental_frequency,
                                        run_output_dir + "/debug_integrand.dat");
      Debug::export_radiation_phase(electron_beam[0], detector->get_point(0), sim_par.fundamental_frequency,
                                    run_output_dir + "/debug_exponent.dat");
    }

    if (IoUtils::get_required(simulation_config, "dense_frequency_spectrum") == "true" &&
        detector->get_total_points() != 1) {
      std::cerr << "Error: dense_frequency_spectrum=true requires a detector collapsed to a single point, "
                << "but this detector has " << detector->get_total_points() << " points -- only the first "
                << "detector point will be used.\n";
      detector->restrict_to_first_point();
    }

    size_t num_threads = IoUtils::get_num_threads(simulation_config);
    double total_cpu_seconds = 0.0;
    auto simulation_start = std::chrono::steady_clock::now();
    Simulation::RadiationField radiation_field = Simulation::run_simulation(
        simulation_config, *laser, *detector, electron_beam, sim_par.frequencies, num_threads, total_cpu_seconds);
    auto simulation_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> simulation_elapsed = simulation_end - simulation_start;
    std::cout << "Simulation time:      " << simulation_elapsed.count() << " s\n";
    std::cout << "Number of particles:  " << electron_beam.size() << "\n";
    std::cout << "Number of screen pts: " << detector->get_total_points() << "\n";
    std::cout << "Number of threads:    " << num_threads << "\n";
    // Simulation::run_simulation itself already printed a per-thread timing breakdown, plus the
    // aggregate "Time per electron"/"Time per screen pt" (summed across every thread's own elapsed
    // time, not derived from this wall-clock simulation_elapsed -- see run_simulation's doc comment
    // on why wall-clock time alone would understate the true per-unit cost by roughly num_threads).

    Radiation::plot_radiation_field(radiation_field, sim_par.frequencies, sim_par.fundamental_frequency,
                                    run_output_dir + "/radiation_field.dat");

    Logging::write_run_log(simulation_config, *laser, *detector, electron_beam.size(), sim_par, num_threads,
                           simulation_elapsed.count(), total_cpu_seconds, run_output_dir + "/run_log.txt");
    std::cout << "Successfully exported run log to " << run_output_dir << "/run_log.txt\n";

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}