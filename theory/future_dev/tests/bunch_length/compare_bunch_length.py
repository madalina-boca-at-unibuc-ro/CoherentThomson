"""Point-by-point comparison of the bunch-length runs (README.md here, "Analysis").

For each sigma > 0 run and harmonic N, against the near-zero-sigma reference on the same transverse sample:
  c = <F_ref, F_sig> / ||F_ref||^2          expect |c| = rho(q_N) = exp(-(q_N sigma)^2 / 2), arg c ~ 0
  mismatch = 1 - |<F_ref, F_sig>|^2 / (||F_ref||^2 ||F_sig||^2)
  residual r = ||F_sig - c F_ref||^2 / ||F_ref||^2, and eta = r / (1 - rho^2): the reference's incoherent fraction
  (should be the same for every run if the residual is the incoherent remainder N<|e|^2>(1 - rho^2)).
F is the total (LR + SR + BR) Faraday tensor, the six independent components; the circular screen is equal-area,
so every point has the same weight.

Usage: python3 compare_bunch_length.py [index.txt [runs_dir]]
  defaults: ~/output/progress_report_20261029/logs/bunch_length/index.txt and .../progress_report_20261029/bunch_length;
  the index records the runs' original paths (~/output_tests/bunch_length/<run>), so a run that has moved is looked up
  by its folder name in runs_dir.
"""

import re
import sys
from pathlib import Path

import numpy as np

REPORT = Path.home() / "output/progress_report_20261029"
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def read_runs(index_path, runs_dir):
  runs = {}
  for line in Path(index_path).read_text().splitlines():
    m = re.search(r"END\s+(\S+)\s+rc=0\s+run=(\S+)", line)
    if m:
      run = Path(m.group(2))
      runs[m.group(1)] = run if run.is_dir() else Path(runs_dir) / run.name
  return runs


def read_sigma_lambda(run):
  for line in (run / "config.cfg").read_text().splitlines():
    tok = line.split()
    if tok and tok[0] == "beam_sigma_z":
      assert tok[2] == "lambda", tok
      return float(tok[1])
  raise KeyError("beam_sigma_z")


def read_field(run):
  path = run / "radiation_field.dat"
  with open(path) as f:
    lines = (line for line in f)
    n_skip = 0
    for line in lines:
      n_skip += 1
      if not line.startswith("#"):
        header = line.split()
        break
  data = np.loadtxt(path, skiprows=n_skip)
  col = {name: i for i, name in enumerate(header)}
  n_omega = int(data[:, col["i_omega"]].max()) + 1
  n_screen = int(data[:, col["i_screen"]].max()) + 1
  omega = data[:: n_screen, col["omega"]]
  F = np.zeros((n_omega, n_screen, len(PAIRS)), dtype=complex)
  for k, (mu, nu) in enumerate(PAIRS):
    for part in ("LR", "SR", "BR"):
      re_ = data[:, col[f"{part}_F{mu}{nu}_re"]]
      im_ = data[:, col[f"{part}_F{mu}{nu}_im"]]
      F[:, :, k] += (re_ + 1j * im_).reshape(n_omega, n_screen)
  return omega, F


def read_k1_over_kL(run):
  # q_N = N k_1 + k_L head-on, on axis; k_s at N = 1 in units of k_L is in the form-factor section of run_log.txt
  m = re.search(r"k_s\s+\(.*?\)\s+a\.u\.\s+\[\(\S+,\s*\S+,\s*(\S+)\)\s*k_L\]", (run / "run_log.txt").read_text())
  return float(m.group(1))


def main():
  index = sys.argv[1] if len(sys.argv) > 1 else REPORT / "logs/bunch_length/index.txt"
  runs_dir = sys.argv[2] if len(sys.argv) > 2 else REPORT / "bunch_length"
  runs = read_runs(index, runs_dir)
  ref = runs.pop("ref")
  xy_ref = np.loadtxt(ref / "electron_beam_scatter.dat", skiprows=2)[:, :2]
  omega, F_ref = read_field(ref)
  k1 = read_k1_over_kL(ref)  # k_1 / k_L
  two_pi = 2 * np.pi
  print(f"reference {ref.name}, k_1/k_L = {k1:.6f}, harmonics {omega}")
  print(f"{'run':6} {'sigma[lam]':>10} {'N':>2} {'q sigma':>8} {'rho':>8} {'|c|':>8} {'|c|/rho':>8} "
        f"{'arg c':>9} {'mismatch':>9} {'resid':>9} {'eta':>9}")
  for name, run in runs.items():
    xy = np.loadtxt(run / "electron_beam_scatter.dat", skiprows=2)[:, :2]
    same_xy = np.array_equal(xy, xy_ref)
    sigma = read_sigma_lambda(run)
    om, F = read_field(run)
    assert np.allclose(om, omega)
    for i, N in enumerate(omega):
      a, b = F_ref[i].ravel(), F[i].ravel()
      qs = (N * k1 + 1.0) * two_pi * sigma  # k_L = 2 pi / lambda
      rho = np.exp(-qs**2 / 2)
      nref = np.vdot(a, a).real
      c = np.vdot(a, b) / nref
      mismatch = 1 - abs(np.vdot(a, b)) ** 2 / (nref * np.vdot(b, b).real)
      resid = np.linalg.norm(b - c * a) ** 2 / nref
      eta = resid / (1 - rho**2)
      print(f"{name:6} {sigma:10.4e} {N:2.0f} {qs:8.4f} {rho:8.5f} {abs(c):8.5f} {abs(c) / rho:8.5f} "
            f"{np.angle(c):9.2e} {mismatch:9.2e} {resid:9.2e} {eta:9.2e}" + ("" if same_xy else "  (x,y DIFFER)"))


if __name__ == "__main__":
  main()
