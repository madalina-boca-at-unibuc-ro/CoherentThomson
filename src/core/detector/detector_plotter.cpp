#include "../include/detector/detector_plotter.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "../include/io_utils/io_utils.hpp"

namespace Core::Detector {
void plot_detector(const Detector_2D& detector, const std::string& type_of_plot, const std::string& output_dir,
                   const double& axes_scale, const std::string& axes_scale_string) {
  if (detector.get_type_name() == "SphericalDetector") {
    const auto* sph_detector = dynamic_cast<const SphericalDetector*>(&detector);
    if (type_of_plot == "stereographic_plot") {
      std::filesystem::path filepath = std::filesystem::path(output_dir) / "detector_stereographic.dat";
      std::ofstream file(filepath);
      if (file.is_open()) {
        file << std::scientific << std::setprecision(6);
        file << "# axes_unit " + axes_scale_string + "\n";
        file << "x_proj y_proj\n";
        size_t total_points = sph_detector->get_total_points();
        for (size_t idx = 0; idx < total_points; ++idx) {
          auto [x_proj, y_proj] = sph_detector->get_stereographic_projection(idx);
          file << x_proj / axes_scale << " " << y_proj / axes_scale << "\n";
        }
        std::cout << "Successfully exported stereographic projection to " << filepath.string() << "\n";
      } else {
        std::cerr << "Failed to open " << filepath.string() << " for writing.\n";
      }
    } else {
      throw std::invalid_argument("plot_detector: unrecognized type_of_plot '" + type_of_plot +
                                  "' for SphericalDetector");
    }
  }
}

void plot_detector_scatter(const Detector_2D& detector, const std::string& output_dir, const double& axes_scale,
                           const std::string& axes_scale_string) {
  std::filesystem::path filepath = std::filesystem::path(output_dir) / "detector_scatter.dat";
  detector.export_coordinates(filepath.string(), axes_scale, axes_scale_string);
  std::cout << "Successfully exported detector scatter coordinates to " << filepath.string() << "\n";
}

}  // namespace Core::Detector