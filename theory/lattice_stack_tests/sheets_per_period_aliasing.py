#!/usr/bin/env python3
"""Minimum sheets per modulation period S before aliasing spoils the moved-sheet sum b_n (no solver runs).

One period of S sheets at phi0_j = (j + 1/2) 2 pi / S, moved to phi0 + zeta sin(phi0), gives the S-point sum of
exp(i n (phi0 + zeta sin phi0 - phase)). Its exact value is (-1)^n J_n(n zeta) e^{-i n phase}; the S-point sum adds the
aliases J_{mS-n}(n zeta), m != 0. The table is the smallest S whose error stays below tol * max(|b_n|, 1e-3) (worst of
phase 0, pi/2, pi; also for the next five S, as the error isn't monotonic), against the rule n (1 + zeta) + 5 n^(1/3)
that Particle::modulation_warnings uses. n = N Lambda_mod / Lambda_1. The floor 1e-3 (of full bunching) keeps the
tolerance meaningful for harmonics the modulation barely excites, where a relative error would demand aliases below
double-precision roundoff.

Run: python3 theory/lattice_stack_tests/sheets_per_period_aliasing.py   (a few seconds; results in
sheets_per_period_aliasing_results.txt)
"""
import numpy as np
from scipy.special import jv

PHASES = (0.0, 0.5 * np.pi, np.pi)


def sheet_sum(n, zeta, sheets, phase):
    phi0 = (np.arange(sheets) + 0.5) * 2 * np.pi / sheets
    return np.mean(np.exp(1j * n * (phi0 + zeta * np.sin(phi0) - phase)))


def minimum_sheets(n, zeta, tol):
    b = (-1) ** n * jv(n, n * zeta)
    exact = [b * np.exp(-1j * n * phase) for phase in PHASES]
    limit = tol * max(abs(b), 1e-3)

    def good(sheets):
        return all(abs(sheet_sum(n, zeta, sheets, p) - e) <= limit for p, e in zip(PHASES, exact))

    return next(s for s in range(1, 5000) if all(good(s2) for s2 in range(s, s + 6)))


print("Minimum sheets per period S for an error of b_n below tol * max(|b_n|, 1e-3); rule = n (1 + zeta) + 5 n^(1/3)")
print("   n   zeta     |b_n|    1e-2  1e-4  1e-6   rule")
for zeta in (0.5, 1.0, 1.111, 1.5):
    for n in (1, 2, 3, 5, 10, 20, 40, 60):
        rule = n * (1 + zeta) + 5 * n ** (1 / 3)
        found = [minimum_sheets(n, zeta, tol) for tol in (1e-2, 1e-4, 1e-6)]
        print(f"{n:4d} {zeta:6.3f}  {abs(jv(n, n * zeta)):.2e}  {found[0]:4d}  {found[1]:4d}  {found[2]:4d}  {rule:6.1f}")
