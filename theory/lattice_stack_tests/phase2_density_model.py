#!/usr/bin/env python3
"""Phase 2: the cold-beam density model, as a standalone reference. It does not touch the solver.

The drift map is phi = phi_0 + zeta sin(phi_0), with phi_0 uniform, so b_n = <exp(-i n phi)> = (-1)^n J_n(n zeta).
Two ways to put it on sheets, both at fixed total charge (sum of weights = 1):

  charges  equidistant sheets, each with the charge integrated over its bin (the brief's route)
  moved    equal charges, each sheet displaced from its uniform position by the map

Geometry: a stack of P periods with S sheets per period has uniform positions z_j = (j - (S P - 1)/2) d, d = Lambda_mod/S,
and phi = K z + origin, K = 2 pi / Lambda_mod. All lengths below are K z, so Lambda_mod is 2 pi.
The solver's convention is b(q) = sum_j w_j exp(+i q z_j), so at q = n K it gives conj(b_n) exp(-i n origin).
For harmonic N of a period ratio r = Lambda_mod / Lambda_1, q_N = N q_1 means n = N r.

From the repository root:
  .venv/bin/python theory/lattice_stack_tests/phase2_density_model.py                 # checkpoint; writes phase2_results.txt
  .venv/bin/python theory/lattice_stack_tests/phase2_density_model.py --zeta 1 --N 2 --ratio 20 --sheets 120
  .venv/bin/python theory/lattice_stack_tests/phase2_density_model.py --zeta 1 --ratio 20 --sheets 280 --q-over-K 20.5
  .venv/bin/python theory/lattice_stack_tests/phase2_density_model.py --run <run_folder>    # b from a run's beam file
"""
import argparse
import re
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.special import jv

from lattice_stack_common import form_factor_b

HERE = Path(__file__).resolve().parent
TWO_PI = 2.0 * np.pi
SITES_PER_SHEET = 1076  # run A's lattice, only used to quote trajectory counts


def bessel_b(n, zeta):
    """The benchmark b_n = (-1)^n J_n(n zeta), for n >= 1."""
    return (-1.0) ** n * jv(n, n * zeta)


def direct_b(n, zeta, points=2_000_000):
    """b_n from the phase integral itself (midpoint rule on a periodic integrand)."""
    phi_0 = (np.arange(points) + 0.5) * TWO_PI / points
    return np.mean(np.exp(-1j * n * (phi_0 + zeta * np.sin(phi_0))))


@lru_cache(maxsize=None)
def exact_bin_charges(zeta, sheets, first_centre=None):
    """Charge per bin of width 2 pi / sheets, from the map's exact crossings of the bin edges (all branches).

    The bins are centred on first_centre + i * 2 pi / sheets; the default is the midpoints (i + 1/2) 2 pi / sheets.
    Between its turning points the map is monotonic, so each branch is cut at the edges it crosses and every piece goes to
    the bin it lands in. Bisection brings the crossings to rounding error, which leaves no quadrature to converge.
    """
    width = TWO_PI / sheets
    if first_centre is None:
        first_centre = width / 2
    first_edge = first_centre - width / 2
    mapped = lambda phi_0: phi_0 + zeta * np.sin(phi_0) - first_edge

    branch_ends = [0.0, TWO_PI]
    if zeta >= 1.0:
        turn = np.arccos(-1.0 / zeta)
        branch_ends = [0.0, turn, TWO_PI - turn, TWO_PI]

    charges = np.zeros(sheets)
    for a, b in zip(branch_ends[:-1], branch_ends[1:]):
        if b <= a:
            continue
        g_a, g_b = mapped(a), mapped(b)
        k = np.arange(np.ceil(min(g_a, g_b) / width), np.floor(max(g_a, g_b) / width) + 1)
        targets = k * width
        low, high = np.full(k.shape, a), np.full(k.shape, b)
        rising = g_b > g_a
        for _ in range(70):
            mid = 0.5 * (low + high)
            go_up = (mapped(mid) < targets) == rising
            low, high = np.where(go_up, mid, low), np.where(go_up, high, mid)
        cuts = np.concatenate(([a], np.sort(0.5 * (low + high)), [b]))
        centres = 0.5 * (cuts[:-1] + cuts[1:])
        bins = np.floor(mapped(centres) / width).astype(np.int64) % sheets
        np.add.at(charges, bins, np.diff(cuts) / TWO_PI)
    return charges


def quadrature_bin_charges(zeta, sheets, points, first_centre=None):
    """The same charges by mapping `points` uniform phases into the bins, the route the brief suggests."""
    width = TWO_PI / sheets
    if first_centre is None:
        first_centre = width / 2
    first_edge = first_centre - width / 2
    counts = np.zeros(sheets)
    chunk = 2_000_000
    for start in range(0, points, chunk):
        phi_0 = (np.arange(start, min(start + chunk, points)) + 0.5) * TWO_PI / points
        bins = np.floor((phi_0 + zeta * np.sin(phi_0) - first_edge) / width).astype(np.int64) % sheets
        counts += np.bincount(bins, minlength=sheets)
    return counts / points


def uniform_positions(sheets, periods):
    """K z_j of the uniform stack, centred on z = 0 as in electron_factory.cpp."""
    total = sheets * periods
    return (np.arange(total) - (total - 1) / 2) * TWO_PI / sheets


def charges_stack(zeta, sheets, periods=1, origin=0.0):
    """Returns (K z_j, w_j) for equidistant sheets carrying the exact bin charges."""
    kz = uniform_positions(sheets, periods)
    first = float(np.mod(kz[0] + origin, TWO_PI))
    return kz, np.tile(exact_bin_charges(zeta, sheets, first), periods) / periods


def moved_stack(zeta, sheets, periods=1, origin=0.0):
    """Returns (K z_j, w_j) for equal charges displaced by the map. Sheets can leave the envelope by up to zeta / K."""
    kz = uniform_positions(sheets, periods)
    return kz + zeta * np.sin(kz + origin), np.full(kz.shape, 1.0 / kz.size)


def nearest_sheet_stack(zeta, sheets, origin, points=40_000_000):
    """Charges of a one-period stack by another route: each uniform phase goes to the sheet nearest its position."""
    kz = uniform_positions(sheets, 1)
    width = TWO_PI / sheets
    counts = np.zeros(sheets)
    chunk = 2_000_000
    for start in range(0, points, chunk):
        phi_0 = (np.arange(start, min(start + chunk, points)) + 0.5) * TWO_PI / points + origin
        position = phi_0 + zeta * np.sin(phi_0) - origin
        counts += np.bincount(np.rint((position - kz[0]) / width).astype(np.int64) % sheets, minlength=sheets)
    return kz, counts / points


def solver_b(kz, weights, q_over_K):
    """b(q) in the solver's convention."""
    return np.sum(weights * np.exp(1j * q_over_K * kz))


def paper_b(zeta, sheets, n, method):
    """The sheet sum over one period, in the paper's convention: <exp(-i n phi)>."""
    kz, weights = method(zeta, sheets)
    return np.sum(weights * np.exp(-1j * n * kz))


def relative(error, reference):
    """Error relative to the benchmark; the absolute error where the benchmark vanishes (zeta = 0)."""
    return abs(error) / abs(reference) if abs(reference) > 1e-12 else abs(error)


def density_at(phi, zeta, points=200_000):
    """The brief's density rho/rho_mean, summing 1/|1 + zeta cos(phi_0)| over the roots of the map at phi."""
    phi_0 = np.linspace(0.0, TWO_PI, points + 1)
    g = phi_0 + zeta * np.sin(phi_0)
    total = 0.0
    for m in range(-3, 4):
        d = g - (phi + m * TWO_PI)
        hits = np.nonzero(d[:-1] * d[1:] < 0)[0]
        for i in hits:
            root = phi_0[i] - d[i] * (phi_0[i + 1] - phi_0[i]) / (d[i + 1] - d[i])
            total += 1.0 / abs(1.0 + zeta * np.cos(root))
    return total


# ---------------------------------------------------------------------------- checkpoint sections


def check_charge(out):
    out("=== 1. Charge conservation and sign, per period (exact crossings; 4e7-point quadrature)")
    out("The exact charges sum to 1 by construction (a partition of phi_0); the quadrature's count every point once.")
    out("min w and max w are relative to the mean charge. Bins here have their edges at multiples of 2 pi / S.")
    out(" zeta  sheets  |sum-1| exact  |sum-1| quad    min w   max w  negative?  trajectories")
    for zeta in (0.0, 1.0, 1.111, 1.5):
        for sheets in (14, 120, 280, 560):
            w = exact_bin_charges(zeta, sheets)
            quad = quadrature_bin_charges(zeta, sheets, 40_000_000)
            out(f"{zeta:5.3f} {sheets:7d}  {abs(w.sum() - 1):13.1e}  {abs(quad.sum() - 1):12.1e}  {w.min() * sheets:7.3f} "
                f"{w.max() * sheets:7.2f}  {str(bool((w < 0).any())):9s}  {sheets * SITES_PER_SHEET}")


def check_density(out):
    out("\n=== 2. Bin charges against the brief's density formula at the bin centre (rho/mean; away from folds)")
    out(" zeta  sheets  median |S w - rho|/rho")
    for zeta in (0.5, 1.5):
        for sheets in (280, 1120, 4480):
            w = exact_bin_charges(zeta, sheets)
            picks = np.arange(0, sheets, sheets // 40)
            centres = (picks + 0.5) * TWO_PI / sheets
            rho = np.array([density_at(c, zeta) for c in centres])
            out(f"{zeta:5.3f} {sheets:7d}  {np.median(np.abs(sheets * w[picks] - rho) / rho):.2e}")


def check_integral(out):
    out("\n=== 3. Direct phase integral against the Bessel benchmark")
    out("  n  zeta  (-1)^n J_n(n zeta)  |direct - Bessel|")
    for n in (1, 2, 20, 40):
        for zeta in (0.0, 1.0, 1.111, 1.5):
            ref = bessel_b(n, zeta)
            out(f"{n:3d} {zeta:5.3f}  {ref:+18.10f}  {abs(direct_b(n, zeta) - ref):.1e}")


def check_sheet_error(out):
    out("\n=== 4. Sheet sum against the benchmark. Errors relative to |b_n| (absolute at zeta = 0)")
    out("charges: spacing = exact bin charges on the sheets, the error of the sheets themselves; after sinc = the same with")
    out("         the bin factor sin(x)/x, x = n pi / S, divided out; integration = 4e6 and 4e7-point quadrature bins")
    out("         against exact bins, the separate error of computing the charges. moved: the sheet sum.")
    out("  n  zeta  sheets  benchmark   | charges: spacing  after sinc   integr. 4e6  integr. 4e7 | moved")
    cases = [(1, [7, 14, 28, 56]), (2, [14, 28, 56]), (20, [60, 80, 120, 280, 560, 1120]), (40, [120, 160, 280, 560, 1120])]
    for n, sheet_list in cases:
        for zeta in (0.0, 1.0, 1.111, 1.5):
            ref = bessel_b(n, zeta)
            for sheets in sheet_list:
                b = paper_b(zeta, sheets, n, charges_stack)
                x = n * np.pi / sheets
                centres = (np.arange(sheets) + 0.5) * TWO_PI / sheets
                sinc_fixed = b / (np.sin(x) / x)
                integration = []
                for points in (4_000_000, 40_000_000):
                    quad = quadrature_bin_charges(zeta, sheets, points)
                    integration.append(abs(np.sum((quad - exact_bin_charges(zeta, sheets)) * np.exp(-1j * n * centres))))
                moved = paper_b(zeta, sheets, n, moved_stack)
                out(f"{n:3d} {zeta:5.3f} {sheets:7d}  {ref:+.6f}  | {relative(b - ref, ref):16.2e} {relative(sinc_fixed - ref, ref):11.2e}"
                    f"  {relative(integration[0], ref):11.1e}  {relative(integration[1], ref):11.1e} | {relative(moved - ref, ref):9.2e}")


def check_other_harmonics(out):
    out("\n=== 5. One stack, every harmonic: error at n = N r for r = 20. Which N does a sheet count still resolve?")
    out(" zeta  sheets  method   " + "".join(f"  N={n // 20:<7d}" for n in range(20, 141, 20)))
    for zeta in (1.0, 1.111, 1.5):
        for sheets in (60, 120, 180, 280, 560):
            for name, method in (("charges", charges_stack), ("moved", moved_stack)):
                row = "".join(f"  {relative(paper_b(zeta, sheets, n, method) - bessel_b(n, zeta), bessel_b(n, zeta)):9.1e}"
                              for n in range(20, 141, 20))
                out(f"{zeta:5.3f} {sheets:7d}  {name:7s} {row}")


def check_convention(out):
    out("\n=== 6. Solver's convention and phase origin: b(nK) against conj(b_n) exp(-i n origin), whole periods")
    out("  n  zeta  sheets  periods  origin  |charges - pred|/|b_n| (exact bins include the bin error)   |moved - pred|/|b_n|")
    for n, zeta, sheets in ((1, 1.0, 56), (20, 1.0, 280), (20, 1.5, 560)):
        for periods in (1, 2, 3):
            for origin in (0.0, 0.7, -2.3):
                pred = np.conj(bessel_b(n, zeta)) * np.exp(-1j * n * origin)
                b_c = solver_b(*charges_stack(zeta, sheets, periods, origin), n)
                b_m = solver_b(*moved_stack(zeta, sheets, periods, origin), n)
                out(f"{n:3d} {zeta:5.3f} {sheets:7d} {periods:8d} {origin:7.2f}  {abs(b_c - pred) / abs(pred):.2e}"
                    f"{'':50s}{abs(b_m - pred) / abs(pred):.2e}")
    out("charges by a separate route, one period: each phase to its nearest sheet (4e7 points) against exact crossings:")
    for n, zeta, sheets, origin in ((20, 1.0, 280, 0.0), (20, 1.0, 280, 0.7), (20, 1.0, 280, -2.3), (20, 1.5, 560, 0.7)):
        b_e = solver_b(*charges_stack(zeta, sheets, 1, origin), n)
        b_q = solver_b(*nearest_sheet_stack(zeta, sheets, origin), n)
        out(f"  n={n} zeta={zeta} S={sheets} origin={origin:+.1f}: exact {b_e:+.6f}, nearest {b_q:+.6f}, difference {abs(b_e - b_q):.1e}")
    out("uniform control (zeta = 0, whole periods), |b(nK)|:")
    for n, sheets, periods in ((1, 14, 1), (20, 280, 1), (20, 280, 3), (40, 280, 1)):
        out(f"  n={n} S={sheets} P={periods}: {abs(solver_b(*charges_stack(0.0, sheets, periods), n)):.1e}")


def fold_into_envelope(stack, periods):
    """Moved sheets folded back into H: the same periodic density, cut at the envelope like the charges' one."""
    half = np.pi * periods
    return np.mod(stack[0] + half, 2 * half) - half, stack[1]


def check_finite_train(out):
    out("\n=== 7. Finite train: |b(q)| around the peak and at the sheet replica, n = 20, P = 3, q = n K + offset K / P")
    n, periods = 20, 3
    for zeta, origin in ((1.0, 0.0), (1.0, 0.7), (1.5, 0.0), (1.5, 0.7)):
        moved_fine = moved_stack(zeta, 5600, periods, origin)
        stacks = {"charges 280": charges_stack(zeta, 280, periods, origin),
                  "charges 5600": charges_stack(zeta, 5600, periods, origin),
                  "charges 56000": charges_stack(zeta, 56000, periods, origin),
                  "moved 120": moved_stack(zeta, 120, periods, origin),
                  "moved 280": moved_stack(zeta, 280, periods, origin),
                  "moved 5600": moved_fine,
                  "folded 5600": fold_into_envelope(moved_fine, periods)}
        spill = max(0.0, np.abs(moved_fine[0]).max() - np.pi * periods) / TWO_PI
        out(f"zeta = {zeta}, origin = {origin}: moved sheets leave H by {spill:.4f} Lambda_mod at each end")
        out("    offset  " + "".join(f"{name:>15s}" for name in stacks))
        for offset in (0.0, 0.5, 1.0, 1.5, 2.0):
            q = n + offset / periods
            out(f"{offset:10.1f}  " + "".join(f"{abs(solver_b(*stack, q)):15.6f}" for stack in stacks.values()))
    out("At q d = 2 pi (q = S K), zeta = 1, |b|: the replica of equidistant sheets, which is 1 for any charges")
    for name, method, sheets in (("charges", charges_stack, 280), ("moved", moved_stack, 280),
                                 ("charges", charges_stack, 120), ("moved", moved_stack, 120)):
        out(f"  {name:7s} S = {sheets}: {abs(solver_b(*method(1.0, sheets, periods), sheets)):.6f}")


def b_from_run(run):
    """Recomputes b from a run's electron_beam_scatter.dat (x y z in lambda, equal weights) and the log's q."""
    data = np.loadtxt(Path(run) / "electron_beam_scatter.dat", skiprows=2)
    text = (Path(run) / "run_log.txt").read_text()
    match = re.search(r"^\s+q\s+\(.*?\[\((\S+), (\S+), (\S+)\) k_L\]", text, re.M)
    if not match:
        raise RuntimeError(f"no q in {run}/run_log.txt")
    q_kL = np.array([float(v) for v in match.groups()])
    return np.mean(np.exp(1j * TWO_PI * (data @ q_kL))), form_factor_b(run), data.shape[0]


def check_runs(out):
    out("\n=== 8. Against the solver's log: b from the run's beam file and the log's 7-digit q (uniform stacks only so far)")
    home = Path.home()
    runs = sorted((home / "output" / "lattice_stack").glob("2*")) + sorted(
        (home / "output_tests" / "lattice_stack" / "model_config_check_scatter_M14").glob("2*"))
    for run in runs:
        if not (run / "electron_beam_scatter.dat").exists() or not (run / "run_log.txt").exists():
            continue
        mine, log, count = b_from_run(run)
        out(f"{run.name}: {count} electrons, tool {mine.real:+.7f}{mine.imag:+.1e}j, log {log.real:+.7f}{log.imag:+.1e}j, "
            f"difference {abs(mine - log):.1e}")


def checkpoint():
    path = HERE / "phase2_results.txt"
    with open(path, "w") as f:
        def out(line):
            print(line, flush=True)
            f.write(line + "\n")
        for section in (check_charge, check_density, check_integral, check_sheet_error, check_other_harmonics,
                        check_convention, check_finite_train, check_runs):
            section(out)
    print(f"\nWritten to {path}")


def single_case(args):
    n = args.N * args.ratio
    q = args.q_over_K if args.q_over_K is not None else float(n)
    print(f"zeta = {args.zeta}, S = {args.sheets} sheets per period, P = {args.periods}, origin = {args.origin}, "
          f"q/K = {q}, {args.sheets * args.periods} sheets, {args.sheets * args.periods * SITES_PER_SHEET} trajectories")
    if args.q_over_K is None and args.zeta >= 0:
        pred = np.conj(bessel_b(n, args.zeta)) * np.exp(-1j * n * args.origin)
        print(f"benchmark  conj(b_{n}) exp(-i n origin) = {pred:+.8f}")
    for name, method in (("charges", charges_stack), ("moved", moved_stack)):
        b = solver_b(*method(args.zeta, args.sheets, args.periods, args.origin), q)
        print(f"{name:8s} b(q) = {b:+.8f}, |b| = {abs(b):.8f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--zeta", type=float, help="drift ratio; with it, compute one case instead of the checkpoint")
    parser.add_argument("--N", type=int, default=1, help="scattered harmonic")
    parser.add_argument("--ratio", type=int, default=1, help="Lambda_mod / Lambda_1")
    parser.add_argument("--sheets", type=int, default=280, help="sheets per modulation period")
    parser.add_argument("--periods", type=int, default=1)
    parser.add_argument("--origin", type=float, default=0.0, help="phase origin, phi = K z + origin")
    parser.add_argument("--q-over-K", type=float, help="any q in units of K, instead of q = N ratio K")
    parser.add_argument("--run", help="a run folder: compare b from its beam file with its log")
    args = parser.parse_args()
    if args.run:
        mine, log, count = b_from_run(args.run)
        print(f"{count} electrons: tool {mine:+.8f}, log {log:+.8f}, difference {abs(mine - log):.1e}")
    elif args.zeta is not None:
        single_case(args)
    else:
        checkpoint()


if __name__ == "__main__":
    sys.exit(main())
