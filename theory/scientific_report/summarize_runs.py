"""Prints the report's per-run numbers from each run folder's run_log.txt (after plot_observables.py and
lg_decomposition.py have written their sections): the Sec. 3.4 angular momentum per photon and the Sec. 4.3 LG
decomposition of E_-.

Run from anywhere: python3 summarize_runs.py <run_folder> [<run_folder> ...]
"""
import re
import sys
from pathlib import Path


def section(text, title):
    """Lines of the run_log.txt section `title` (up to the next underlined heading)."""
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(title)) + 2
    end = next((i - 1 for i in range(start + 1, len(lines)) if re.fullmatch(r'-{5,}', lines[i])), len(lines))
    return lines[start:end]


def table_rows(lines):
    """Numeric rows (i_omega first) of a whitespace table."""
    return [list(map(float, l.split())) for l in lines if re.match(r'\s+\d+\s+[\d.]+\s', l)]


def summarize(run):
    text = (run / 'run_log.txt').read_text()
    omega_1 = float(re.search(r'fundamental_frequency\s+(\S+) a\.u\.', text).group(1))
    count = re.search(r'particle count\s+(\d+)', text).group(1)
    lattice = re.search(r'polar_lattice.*', text)
    print(f'{run}  ({lattice.group(0).strip() if lattice else f"random, {count} electrons"})')

    # Sec. 3.4: spin and orbital angular momentum per photon, omega_N * flux / P_z (hbar = 1)
    obs = table_rows(section(text, 'Screen-integrated observables'))
    print('  N   spin    orbital  total')
    for row in obs:
        N, sigma_zz, lambda_zz, p_z = row[1], row[5], row[6], row[9]
        spin, orbital = N * omega_1 * sigma_zz / p_z, N * omega_1 * lambda_zz / p_z
        print(f'  {N:.0f}  {spin:+.3f}  {orbital:+.3f}  {spin + orbital:+.3f}')

    # Sec. 4.3: E_- at the (w0, z0) with the largest mode overlap
    lg = '\n'.join(section(text, 'LG decomposition'))
    print('  N   l   purity  best mode  share  w0[lam]  z0[lam]  sum|a|^2  beyond0.9R  series w0, z0')
    for block in re.split(r'\n\s*\n', lg):
        head = re.search(r'N = (\d+)\), component minus', block)
        if not head:
            continue
        best = re.search(r'z0 = (\S+) lambda_L along propagation \(scan.*?\n\s*w0 = (\S+) lambda_L, where \(p,l\) = '
                         r'\((\d+),(-?\d+)\) holds (\S+) of the field.*?sum\|a\|\^2/sum A\|E\|\^2 = (\S+)', block, re.S)
        purity = re.search(r'azimuthal power: l = (-?\d+): (\S+);', block)
        beyond = re.search(r'power beyond 0\.9 x the screen.s inscribed radius = (\S+)', block)
        series = re.search(r'series waist .*?w0 = (\S+) lambda_L, waist at z0 = (\S+) lambda_L', block)
        z0, w0, p, l, share, total = best.groups()
        flag = ' UNDER-RESOLVED' if 'UNDER-RESOLVED' in block else ''
        print(f'  {head.group(1)}  {purity.group(1):>3}  {float(purity.group(2)):.3f}   LG_{p}^{l:<3}    '
              f'{float(share):.3f}  {float(w0):6.1f}  {float(z0):8.0f}  {float(total.rstrip(";")):.2f}      '
              f'{float(beyond.group(1).rstrip(":")):.4f}      {float(series.group(1)):.2f}, {float(series.group(2)):.3g}'
              f'{flag}')
    print()


if __name__ == '__main__':
    for arg in sys.argv[1:]:
        summarize(Path(arg).expanduser())
