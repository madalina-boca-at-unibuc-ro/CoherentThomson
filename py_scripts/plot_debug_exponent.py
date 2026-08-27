import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file, get_laser_period

def plot_debug_exponent(filepath):
    """
    Plots the phase factor exp(i*(x[0]+R)*k), k = omega/c, vs. tau/T (T = 2*pi/omega, the laser
    period, read from the run's own config.cfg), read from debug_exponent.dat
    (Core::Debug::export_radiation_phase, only written when the config's 'debug' key is true) -- the
    one oscillatory factor common to every long-range/short-range bivector component in
    debug_integrand.dat (plot_debug_integrand.py), factored out into its own file since it's
    identical across all 6 components x 2 (long/short).

    Top panel: the raw (unwrapped) phase in radians vs. tau/T -- its local slope is the instantaneous
    frequency seen in the retarded time x[0]+R, useful for checking the trajectory's tau-sampling
    resolves the phase (i.e. doesn't step by more than ~pi per point). Bottom panel: Re/Im
    (cos/sin) of the phase factor itself, overlaid.
    """
    try:
        data = pd.read_csv(filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    required_cols = {'tau', 'phase', 'exp_re', 'exp_im'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain columns {sorted(required_cols)}")
        sys.exit(1)

    T = get_laser_period(os.path.dirname(filepath))
    data = data.sort_values('tau')
    tau_over_T = data['tau'].to_numpy() / T

    fig, (ax_phase, ax_exp) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax_phase.plot(tau_over_T, data['phase'].to_numpy(), color='darkorange')
    ax_phase.set_title(r"Unwrapped phase: $(x^0 + R) \cdot k$, $k = \omega/c$")
    ax_phase.set_ylabel("phase (rad)")

    ax_exp.plot(tau_over_T, data['exp_re'].to_numpy(),
               label=r"$\mathrm{Re}(e^{i\,\mathrm{phase}}) = \cos(\mathrm{phase})$")
    ax_exp.plot(tau_over_T, data['exp_im'].to_numpy(),
               label=r"$\mathrm{Im}(e^{i\,\mathrm{phase}}) = \sin(\mathrm{phase})$")
    ax_exp.set_title(r"Phase factor $e^{i\,\mathrm{phase}}$")
    ax_exp.set_ylabel("amplitude")
    ax_exp.set_xlabel(r"$\tau / T$  ($T = 2\pi/\omega$, laser period)")
    ax_exp.set_ylim(-1.05, 1.05)
    ax_exp.legend()

    fig.suptitle(f"Radiation integrand phase vs. $\\tau/T$ ({len(tau_over_T)} trajectory points)", fontsize=13,
                fontweight='bold')
    plt.tight_layout()

    png_dir = os.path.join(os.path.dirname(filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)
    output_img = os.path.join(png_dir, "debug_exponent.png")
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        input_file = sys.argv[1]
    else:
        try:
            input_file = find_latest_output_file("debug_exponent.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_debug_exponent.dat>")
            sys.exit(1)
        print(f"No file given; using latest run's debug exponent: {input_file}")

    plot_debug_exponent(input_file)
