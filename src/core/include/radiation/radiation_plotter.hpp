#pragma once

#include <string>
#include <vector>

#include "../simulation/simulation.hpp"

namespace Core::Radiation {

// Exports every element of the coherently-summed long_range/short_range/boundary Faraday tensors
// (boundary is identically zero for radiation_formula="direct" -- see Radiation::compute_radiation's
// doc comment), one row per (frequency, screen-point) pair, for later inspection/plotting. The
// exported "omega" column is
// frequencies_list normalized by fundamental_frequency (Simulation::simulation_parameters::fundamental_frequency,
// i.e. PhysUtils::non_linear_Thomson_formula(k1, p, n2, 1)) rather than left in raw atomic-unit omega, so a
// value of 3.0 reads as "third harmonic" regardless of the observation direction's Doppler shift.
void plot_radiation_field(const Simulation::RadiationField& field, const std::vector<double>& frequencies_list,
                          double fundamental_frequency, const std::string& filepath);

// Exports the INCIDENT laser beam's own complex ("phasor") Faraday tensor -- not the scattered
// radiation -- in the exact same file format plot_radiation_field writes, so it can be compared
// directly (via py_scripts/radiation/plot_field.py's/plot_spherical_components.py's --incident flag)
// against a run's own radiation_field.dat. This gives a clean, exactly-known reference field (no
// coherent-sum noise, no far-field Fresnel/aliasing issues) to sanity-check the scattered field's
// own numerical behavior against.
//
// Evaluated on the plane z_loc = Z of the laser's own frame (epsilon_1, epsilon_2, n), Z = the detector's own
// distance: downstream of the waist, along the laser, whatever the detector's actual orientation. Both
// positions AND components are in that frame, so the file shows the laser unrotated, as if it propagated
// along +Oz (unlike radiation_field.dat's lab-frame components). The grid has the detector's own point counts
// and point order (same i_screen layout as radiation_field.dat), but its size follows the beam
// (incident_screen): for a Laguerre-Gauss laser, half-width 2 * IoUtils::lg_beam_radius_in_w0(p, l) * w(Z), so
// the mode is resolved at any Z; for a plane wave, the detector's own grid. Evaluated at the middle of the
// pulse's flat-top plateau (phi = (get_phi_min()+get_phi_max())/2, i.e. envelope = 0, peak amplitude; ct is
// shifted by Z so the plane sees the same point of the pulse as the waist would) -- see
// LaserField::get_complex_faraday_tensor's own doc comment for why the exact instant chosen doesn't matter.
// Before 2026-10-27 this was the waist plane on the detector's own grid. Only rectangular and circular
// detectors have a well-defined flat local (x, y) plane -- throws for any other type.
void export_incident_field_fourier(const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                                   double fundamental_frequency, const std::string& filepath);

// The plane and grid size export_incident_field_fourier uses. half_width is the circular grid's R_max (R_min =
// 0) or the rectangular grid's x/y half-width (centred), in atomic units; 0 means "the detector's own grid"
// (plane-wave laser). Throws for detectors other than rectangular/circular.
struct IncidentScreen {
  double distance;
  double half_width;
};
IncidentScreen incident_screen(const Laser::LaserField& laser, const Detector::Detector_2D& detector);

// Writes <run_output_dir>/incident_config.cfg: override lines for the incident screen (its distance and size
// in a.u., and detector direction (0, 0) since incident_field.dat is in the laser's frame), followed by a
// verbatim copy of the run's config.cfg. The Python scripts read the first occurrence of a key, so for
// incident_field.dat they read this file instead of config.cfg and need no other special case. Not a config
// for the solver.
void write_incident_config(const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                           const std::string& run_output_dir);

}  // namespace Core::Radiation
