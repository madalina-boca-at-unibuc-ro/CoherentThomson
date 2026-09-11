"""
Shared helper for adding a supplementary Laguerre-Gauss beam-waist (w0) length scale to a
heatmap's axes, alongside its primary axes (already in the run's own length unit, conventionally
'lambda') -- used by every script here that renders a spatial (not angular) length-unit heatmap:
plot_field_heatmap_z0.py, plot_radiation_field.py, plot_spherical_field_components.py,
plot_angular_momentum_flux.py.

Kept as its own leaf module (no import of plot_radiation_field.read_config_value) so every one of
those scripts can import it without a circular-import issue: plot_field_heatmap_z0.py already
imports from plot_radiation_field.py, and plot_radiation_field.py is itself one of this module's
callers -- if this module imported plot_radiation_field.read_config_value at its own top level (or
vice versa), that pair would import each other. _read_config_value below is therefore a small
standalone duplicate of that same 'key value [unit]' line reader, not a re-export.
"""
import os


def _read_config_value(key, config_path):
    """Minimal standalone 'key value [unit]' line reader (see plot_radiation_field.read_config_value
    for the canonical version, used everywhere else -- duplicated, not imported, per the module
    docstring above)."""
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split()
            if parts[0] == key and len(parts) >= 2:
                return parts[1], (parts[2] if len(parts) >= 3 else '')
    raise ValueError(f"No '{key}' key found in '{config_path}'")


def get_laser_lg_w0_in_axes_units(filepath, axes_unit):
    """
    Returns laser_lg_w0 (the Laguerre-Gauss beam waist) expressed in axes_unit (the heatmap's own
    length unit, conventionally 'lambda'), or None if this run's laser_type isn't 'laguerre_gauss'
    (no w0 to label with) or laser_lg_w0's own config unit doesn't match axes_unit (defensive:
    laser_lg_w0 is always given in 'lambda' units by this project's own config convention, matching
    every detector length key's own convention -- but degrades to None instead of silently
    mislabeling the axis if that ever changes, rather than doing a full unit conversion for a case
    that shouldn't occur in practice).

    Reads the run's own config.cfg snapshot next to `filepath`
    (Core::IoUtils::copy_config_to_run_directory), not the live repo config, per the convention
    every other script here already follows.
    """
    config_path = os.path.join(os.path.dirname(filepath), 'config.cfg')
    if not os.path.exists(config_path):
        return None
    try:
        laser_type, _ = _read_config_value('laser_type', config_path)
        if laser_type != 'laguerre_gauss':
            return None
        w0_str, w0_unit = _read_config_value('laser_lg_w0', config_path)
    except ValueError:
        return None
    if w0_unit != axes_unit:
        print(f"Warning: laser_lg_w0's unit ('{w0_unit}') does not match the heatmap's axes_unit "
              f"('{axes_unit}') -- skipping the w0-unit axis labels.")
        return None
    return float(w0_str)


def add_w0_secondary_axes(ax, w0):
    """
    Adds secondary top/right axes in units of the Laguerre-Gauss beam waist w0, alongside the
    primary bottom/left axes (already in the heatmap's own length unit) -- both shown at once
    (rather than replacing one with the other) so either the physical position or its w0-multiple
    can be read directly off the plot. No-op if w0 is None or non-positive (not a laguerre_gauss
    run, a unit mismatch -- see get_laser_lg_w0_in_axes_units -- or an angular, non-length axis,
    which callers should already have excluded before calling this).
    """
    if w0 is None or w0 <= 0:
        return
    secax_x = ax.secondary_xaxis('top', functions=(lambda x: x / w0, lambda x: x * w0))
    secax_x.set_xlabel("$x / w_0$", fontsize=9)
    secax_y = ax.secondary_yaxis('right', functions=(lambda y: y / w0, lambda y: y * w0))
    secax_y.set_ylabel("$y / w_0$", fontsize=9)
