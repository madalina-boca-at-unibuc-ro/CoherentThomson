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
// radiation -- in the exact same file format plot_radiation_field writes, so the same downstream
// angular-momentum analysis (py_scripts/radiation/plot_angular_momentum_*.py) can be pointed at either
// file. This gives a clean, exactly-known reference field (no coherent-sum noise, no far-field
// Fresnel/aliasing issues) to validate those formulas/derivatives against, independent of the
// scattered field's own numerical behavior -- see CLAUDE.md's "analytic incident-field cross-check"
// note for what this was used to find.
//
// Always evaluated at the laser's own canonical-frame beam waist (z_loc=0), regardless of the
// detector's actual configured distance -- and always in the canonical frame (there is no lab-frame
// counterpart, since the detector's actual distance/orientation is intentionally bypassed; unlike
// plot_radiation_field's output this is independent of print_field_in_canonical_frame). Reuses the
// detector's own local 2D grid (get_row_coordinate/get_col_coordinate, ignoring to_lab_frame's
// rotation/distance entirely) so the incident-field screen has identical size/resolution/i_screen
// ordering to the actual radiation_field.dat from the same run -- directly comparable side by side.
// Evaluated at the middle of the pulse's flat-top plateau (phi = (get_phi_min()+get_phi_max())/2,
// i.e. envelope = 0, peak amplitude) -- see LaserField::get_complex_faraday_tensor's own doc comment
// for why the exact instant chosen cannot affect any of the bilinear angular-momentum formulas
// downstream. Only rectangular and circular detectors have a well-defined flat local (x, y) plane to
// reuse (matching the angular-momentum theory docs' own restriction) -- throws for any other type.
void export_incident_field_fourier(const Laser::LaserField& laser, const Detector::Detector_2D& detector,
                                   double fundamental_frequency, const std::string& filepath);

}  // namespace Core::Radiation
