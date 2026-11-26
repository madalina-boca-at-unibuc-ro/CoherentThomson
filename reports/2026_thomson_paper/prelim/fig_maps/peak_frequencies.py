#!/usr/bin/env python3
"""Peak frequencies for the Sec. IV.C field-map runs, from the fig_spectrum line runs.

For each case: the largest maximum of the screen-integrated photon spectrum dN/d(omega/omega_1) (parabolic refinement
over the three grid points around it), and, at that frequency, how the flux d^2N/domega dOmega is spread in theta/theta_*
(fraction beyond 8, 12, 16, 20 theta_*) to size the 2D screen.

Usage: python3 peak_frequencies.py   (reads ~/output/paper_2026/prelim/fig_spectrum/)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'figures'))
sys.argv = sys.argv[:1]
import fig_spectrum as fs  # noqa: E402
import plot_spectrum_map as psm  # noqa: E402

CASES = ['a0_1p0_Nc8_N1', 'a0_1p0_Nc8_N2']


def main():
    print(f'{"run":<24}{"omega_peak":>11}{"grid max":>10}{"beyond 8":>10}{"12":>9}{"16":>9}{"20":>9}  theta 10/50/90%')
    for case in CASES:
        for laser in ('10um', '800nm'):
            r = fs.load(f'{laser}_{fs.CASES[case][1]}')
            w, dN = r['omega'], r['dN']
            k = int(np.nanargmax(dN))
            y0, y1, y2 = dN[k - 1], dN[k], dN[k + 1]
            shift = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            peak = w[k] + shift * (w[k + 1] - w[k])
            t = r['theta'] / r['theta_star']
            flux = r['photons'][k] * psm.annulus_areas(r['theta'])
            cf = np.cumsum(flux) / np.nansum(flux)
            beyond = [1 - cf[np.searchsorted(t, x) - 1] for x in (8, 12, 16, 20)]
            q = [t[np.searchsorted(cf, x)] for x in (0.1, 0.5, 0.9)]
            print(f'{laser + "_" + case:<24}{peak:>11.4f}{w[k]:>10.3f}' + ''.join(f'{b:>9.1e}' for b in beyond)
                  + f'  {q[0]:.1f} / {q[1]:.1f} / {q[2]:.1f}')


if __name__ == '__main__':
    main()
