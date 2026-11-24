#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/include/build_info/build_info.hpp"
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
#include "../core/include/simulation/form_factor.hpp"
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

    // Optional (default false), unlike most keys, so older configs keep working. Checked before the run directory
    // exists, so a refused run leaves nothing behind.
    std::cout << "Build: " << BuildInfo::describe() << "\n";
    auto require_clean = simulation_config.find("require_clean_build");
    if (require_clean != simulation_config.end() && require_clean->second == "true") {
      BuildInfo::require_clean_build();
    }

    std::string run_output_dir = IoUtils::make_run_output_directory(simulation_config);
    std::cout << "Writing run output to " << run_output_dir << "\n";
    IoUtils::copy_config_to_run_directory(config_path, run_output_dir);
    if (BuildInfo::git_dirty) {
      BuildInfo::write_source_patch(run_output_dir + "/source.patch");
      std::cout << "Warning: built from uncommitted changes; wrote them to " << run_output_dir << "/source.patch\n";
    }

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
      // plane (the laser's own transverse plane through the origin) unless the beam's mean momentum
      // is along the laser axis -- centering the window would need reproducing that rotation here,
      // so the window stays beam-radius-sized but always centered on the laser axis.
      // Falls back to the previous windowing (fixed field_heatmap_x/y_min/max for plane_wave,
      // +-2*w0 for laguerre_gauss, whose actual transverse scale is set by its own waist rather than
      // any beam property) only when the beam has no meaningful transverse extent to size off of
      // (radius == 0, e.g. a single on-axis point beam) --
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

      // The instant the plateau start (the peak, with no plateau) crosses z = 0: phi = delay
      // (laser_delay, theory/specs/14_collision_timing.md).
      double t_heatmap = laser->get_delay() / laser->get_omega();

      Laser::export_field_heatmap_z0(*laser, t_heatmap, x_min, x_max, y_min, y_max, heatmap_Nx, heatmap_Ny,
                                     heatmap_axes_scale, heatmap_axes_unit,
                                     run_output_dir + "/laser_field_heatmap_z0.dat");
      std::cout << "Successfully exported laser field heatmap to " << run_output_dir << "/laser_field_heatmap_z0.dat\n";
    }

    Simulation::simulation_parameters sim_par = Simulation::init_simulation_parameters(simulation_config, *laser);

    // Before the beam: a bunched stack's period follows the detector's normal.
    auto detector = Detector::create_detector(simulation_config);

    auto electron_beam = Particle::generate_cylinder_beam(
        simulation_config, sim_par.d_tau_traj, sim_par.simulation_length, *laser, detector->get_normal_direction());

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

    std::string detector_axes_unit = "lambda";
    auto detector_axes_scale = IoUtils::convert_unit_to_number(detector_axes_unit, simulation_config);
    Detector::plot_detector(*detector, "stereographic_plot", run_output_dir, detector_axes_scale, detector_axes_unit);
    if (IoUtils::get_required(simulation_config, "plot_detector_scatter") == "true") {
      Detector::plot_detector_scatter(*detector, run_output_dir, detector_axes_scale, detector_axes_unit);
    }

    // Initial-position form factor at the screen center's Thomson fundamental
    // (theory/specs/06_form_factor_diagnostic.md).
    Simulation::FormFactor form_factor = Simulation::compute_form_factor(electron_beam, *laser, *detector);
    std::cout << "Initial form factor b = (" << form_factor.b.real() << ", " << form_factor.b.imag()
              << "), |b| = " << std::abs(form_factor.b) << "\n";

    // A dense spectrum runs on one point or on a line of points (rectangular Nx or Ny = 1, circular N_phi = 1;
    // theory/specs/15_spectrum_map.md). A 2D grid or a spherical screen keeps only its first point.
    if (IoUtils::get_required(simulation_config, "dense_frequency_spectrum") == "true") {
      if (detector->get_total_points() != 1 && !detector->is_line()) {
        std::cerr << "Error: dense_frequency_spectrum=true needs a single point or a line of points (rectangular "
                  << "Nx or Ny = 1, circular N_phi = 1), but this detector is a " << detector->get_cols() << " x "
                  << detector->get_rows() << " grid -- only the first detector point will be used.\n";
        detector->restrict_to_first_point();
      }
      const Simulation::FieldMemoryEstimate memory =
          Simulation::estimate_field_memory(IoUtils::get_num_threads(simulation_config), electron_beam.size(),
                                            detector->get_total_points(), sim_par.frequencies.size());
      std::cout << "Dense spectrum: " << detector->get_total_points() << " screen points x "
                << sim_par.frequencies.size() << " frequencies; buffers " << memory.total_bytes / 1e9 << " GB ("
                << memory.threads << " threads x " << memory.per_thread_bytes / 1e6 << " MB + "
                << memory.reconstructed_bytes / 1e6 << " MB), radiation_field.dat ~" << memory.file_bytes / 1e6
                << " MB\n";
      if (memory.total_bytes > Simulation::kDenseSpectrumMemoryCap) {
        throw std::runtime_error("dense spectrum too large: " + std::to_string(memory.total_bytes / 1e9) +
                                 " GB of field buffers (" + std::to_string(memory.threads) + " threads x " +
                                 std::to_string(memory.per_thread_bytes / 1e6) + " MB + " +
                                 std::to_string(memory.reconstructed_bytes / 1e6) + " MB) is above the " +
                                 std::to_string(Simulation::kDenseSpectrumMemoryCap / 1e9) +
                                 " GB cap; reduce the screen points, N_omega or num_threads");
      }
    }

    // A first run log before the long part, marked in progress, so a stopped run can still be plotted (beam, laser,
    // detector). Rewritten in full below when the run ends.
    const std::string run_log_path = run_output_dir + "/run_log.txt";
    Logging::write_run_log(simulation_config, *laser, *detector, electron_beam, sim_par, form_factor, std::nullopt,
                           run_log_path);
    std::cout << "Run log (in progress) written to " << run_log_path << "\n";

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

    // Additive diagnostic export: the incident laser beam's own field (not the scattered radiation),
    // in the same file format, as a clean known-analytic reference to compare the scattered field
    // against -- see Radiation::export_incident_field_fourier's doc comment. Only meaningful for
    // rectangular/circular detectors (the only ones with a well-defined flat local (x, y) plane);
    // skipped with a message rather than aborting the run otherwise.
    if (detector->get_type_name() == "RectangularDetector" || detector->get_type_name() == "CircularDetector") {
      Radiation::export_incident_field_fourier(*laser, *detector, sim_par.fundamental_frequency,
                                               run_output_dir + "/incident_field.dat");
      Radiation::write_incident_config(*laser, *detector, run_output_dir);
    } else {
      std::cout << "Skipping incident field export: detector_type '" << detector->get_type_name()
                << "' has no flat local (x, y) plane (only rectangular/circular are supported).\n";
    }

    Logging::write_run_log(simulation_config, *laser, *detector, electron_beam, sim_par, form_factor,
                           Logging::RunTiming{num_threads, simulation_elapsed.count(), total_cpu_seconds}, run_log_path);
    std::cout << "Successfully exported run log to " << run_log_path << "\n";

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}