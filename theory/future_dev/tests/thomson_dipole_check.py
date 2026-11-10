"""Thomson-dipole benchmark check for a thomson_dipole.cfg run (one electron at rest, weak field, spherical screen).

Expected angular pattern at N = 1 (dipole limit): |E|^2 proportional to |A|^2 - |n.A|^2, with the dipole direction
A = zeta_1 epsilon_1 + zeta_2 epsilon_2 read from the run's run_log.txt. For circular polarization along z this is
1 + cos^2(theta); for linear along x, 1 - sin^2(theta) cos^2(phi).

Usage: python3 thomson_dipole_check.py <run_folder>
"""

import re
import sys
from pathlib import Path

import numpy as np

C = 137.035999084  # PhysUtils::AtomicUnits::c


def read_log_vectors(log_text):
    def vec(name):
        m = re.search(rf"^\s*{name}\s+\(([^)]*)\)", log_text, re.M)
        return np.array([float(v) for v in m.group(1).split(",")])

    def cplx(name):
        re_, im_ = vec(name)
        return complex(re_, im_)

    return vec("epsilon_1"), vec("epsilon_2"), cplx("zeta_1"), cplx("zeta_2")


def read_config(run):
    cfg = {}
    for line in (run / "config.cfg").read_text().splitlines():
        line = line.split("#", 1)[0].split()
        if len(line) >= 2:
            cfg[line[0]] = line[1:]
    return cfg


def main():
    run = Path(sys.argv[1]).expanduser()
    cfg = read_config(run)
    assert cfg["detector_type"][0] == "spherical", "needs a spherical detector"
    n_theta = int(cfg["spherical_detector_N_theta"][0])
    n_phi = int(cfg["spherical_detector_N_phi"][0])
    for key in ("spherical_detector_theta_min", "spherical_detector_theta_max", "spherical_detector_phi_min"):
        assert float(cfg[key][0]) in (0.0, 1.0), key
    assert float(cfg["detector_direction_theta"][0]) == 0.0, "expects the screen's axis along +Oz"

    eps1, eps2, z1, z2 = read_log_vectors((run / "run_log.txt").read_text())
    dipole = z1 * eps1 + z2 * eps2

    data = np.loadtxt(run / "radiation_field.dat", comments="#", skiprows=4)
    header = (run / "radiation_field.dat").read_text().splitlines()[3].split()
    col = {name: i for i, name in enumerate(header)}

    def total(mu, nu):
        out = np.zeros(len(data), dtype=complex)
        for rng in ("LR", "SR", "BR"):
            out += data[:, col[f"{rng}_F{mu}{nu}_re"]] + 1j * data[:, col[f"{rng}_F{mu}{nu}_im"]]
        return out

    E = np.stack([C * total(i, 0) for i in (1, 2, 3)], axis=1)
    B = np.stack([total(3, 2), total(1, 3), total(2, 1)], axis=1)

    i_omega = data[:, col["i_omega"]].astype(int)
    i_screen = data[:, col["i_screen"]].astype(int)
    cos_t = 1.0 - 2.0 * (i_screen // n_phi) / (n_theta - 1)
    sin_t = np.sqrt(np.clip(1.0 - cos_t**2, 0.0, None))
    phi = 2.0 * np.pi * (i_screen % n_phi) / (n_phi - 1)
    n = np.stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t], axis=1)

    expected = np.vdot(dipole, dipole).real - np.abs(n @ dipole) ** 2
    print(f"run: {run}")
    print(f"dipole direction A = {np.round(dipole, 6)}  (|A|^2 = {np.vdot(dipole, dipole).real:.6f})")

    peak_n1 = None
    for k in sorted(set(i_omega)):
        sel = i_omega == k
        e2 = np.sum(np.abs(E[sel]) ** 2, axis=1)
        b2 = C**2 * np.sum(np.abs(B[sel]) ** 2, axis=1)
        e_r = np.abs(np.sum(E[sel] * n[sel], axis=1))
        b_r = np.abs(np.sum(B[sel] * n[sel], axis=1))
        omega = data[sel, col["omega"]][0]
        if peak_n1 is None:
            peak_n1 = e2.max()
        print(f"\nharmonic index {k} (omega/omega_1 = {omega:g}): max |E|^2 = {e2.max():.6e}, "
              f"ratio to N = 1 peak = {e2.max() / peak_n1:.3e}")
        if k != 0:
            continue
        ex = expected[sel]
        scale = np.sum(e2 * ex) / np.sum(ex * ex)  # least-squares fit of e2 = scale * expected
        rel = (e2 - scale * ex) / (scale * ex.max())
        print(f"  fit |E|^2 = s * (|A|^2 - |n.A|^2): s = {scale:.6e}")
        print(f"  max |data - fit| / fit peak       = {np.abs(rel).max():.3e}")
        print(f"  c^2|B|^2 vs |E|^2, max rel diff   = {np.max(np.abs(b2 - e2)) / e2.max():.3e}")
        print(f"  |E.n| / |E| (max over screen)     = {np.max(e_r / np.sqrt(e2.max())):.3e}")
        print(f"  |B.n| / |B| (max over screen)     = {np.max(b_r * C / np.sqrt(e2.max())):.3e}")
        print("  theta/pi   cos(theta)   <|E|^2>_phi/s   expected   phi spread/s")
        e2_grid = e2.reshape(n_theta, n_phi)[:, :-1]  # last phi duplicates phi = 0
        ex_grid = ex.reshape(n_theta, n_phi)[:, :-1]
        for i in range(0, n_theta, max(1, (n_theta - 1) // 8)):
            c = 1.0 - 2.0 * i / (n_theta - 1)
            row = e2_grid[i] / scale
            print(f"  {np.arccos(c) / np.pi:7.4f}   {c:+.4f}      {row.mean():.6f}       {ex_grid[i].mean():.6f}   "
                  f"{row.max() - row.min():.2e}")

        absolute_check(run, cfg, E[sel], n[sel])


def absolute_check(run, cfg, E, n):
    """Absolute amplitude against the nonrelativistic far field of the solver's own trajectory (electron.dat):
    E(omega) = -n x (n x A) e^{ikR} / (c^2 R), A = (1/2pi) int dt e^{i omega (t - n.r/c)} d^2r/dt^2
    (atomic units, q = -1, 4 pi epsilon_0 = 1; the theory doc's (1/2pi) int dt e^{i omega t} convention)."""
    log = (run / "run_log.txt").read_text()
    omega = float(re.search(r"^\s*fundamental_frequency\s+(\S+)", log, re.M).group(1))
    lam = 2.0 * np.pi * C / float(cfg["laser_frequency"][0])
    assert cfg["spherical_detector_radius"][1] == "lambda"
    R = float(cfg["spherical_detector_radius"][0]) * lam

    tr = np.loadtxt(run / "electron.dat", skiprows=2)
    tr = tr[tr[:, 0] == 0]
    t, r = tr[:, 2] / C, tr[:, 3:6]
    gamma = tr[:, 6] / C
    accel = tr[:, 11:14] / gamma[:, None] ** 2  # (dp/dtau)/gamma^2 = d^2r/dt^2 to O(beta^2), m = 1
    w = np.full(len(t), t[1] - t[0])
    w[0] = w[-1] = w[0] / 2

    phase = np.exp(1j * omega * (t[None, :] - (n @ r.T) / C))  # (screen, tau)
    A = (phase * w) @ accel / (2.0 * np.pi)
    E_dip = -np.cross(n, np.cross(n, A)) * np.exp(1j * omega / C * R) / (C**2 * R)

    norm_dip = np.linalg.norm(E_dip, axis=1)
    keep = norm_dip > 0.1 * norm_dip.max()  # skip nodes (linear polarization), where the ratio is ill-defined
    ratio = np.linalg.norm(E, axis=1)[keep] / norm_dip[keep]
    dphi = np.angle(np.sum(np.conj(E_dip) * E, axis=1))[keep]
    print(f"  absolute, vs dipole field of electron.dat ({keep.sum()} of {len(keep)} points, nodes skipped):")
    print(f"    |E| / |E_dipole| in [{ratio.min():.6f}, {ratio.max():.6f}]")
    print(f"    phase difference in [{dphi.min():.3e}, {dphi.max():.3e}] rad (a constant is the rounded omega in kR)")


if __name__ == "__main__":
    main()
