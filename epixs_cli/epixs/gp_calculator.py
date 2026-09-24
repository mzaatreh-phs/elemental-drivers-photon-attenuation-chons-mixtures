"""
G-P (Geometric Progression) buildup factor calculator.

Implements the ANSI/ANS-6.4.3-1991 G-P fitting function (Harima's formula):

    B(E, x) = 1 + (b - 1)(K^x - 1) / (K - 1)     for K != 1
    B(E, x) = 1 + (b - 1) x                       for K == 1

    K(E, x) = c * x^a + d * [tanh(x/Xk - 2) - tanh(-2)] / [1 - tanh(-2)]

where x is the penetration depth in mean free paths (mfp) and b is the
buildup factor at 1 mfp. (b, c, a, Xk, d) are looked up per element/mixture
in gp_data.GP_COEFFS and, for an arbitrary compound/mixture, interpolated in
log(Z) between the two standard elements bracketing the material's
equivalent atomic number Zeq(E) -- see zeff.py for Zeq determination.

Validated against EPICS2017's own Lead (Z=82) EBF/EABF output across the
FULL exported grid (23 energies x 18 mfp depths x 2 quantities = 828
points, 30 keV-15 MeV, 1-40 mfp): worst error 0.0005%, median ~0.0001% (see
gp_data.py docstring for the full provenance/validation writeup).
"""

from __future__ import annotations

import math

from gp_data import GP_COEFFS

TANH_M2 = math.tanh(-2.0)


def gp_K(x: float, c: float, a: float, d: float, Xk: float) -> float:
    """K(E, x) -- the G-P formula's shape parameter at depth x mfp."""
    return c * x**a + d * (math.tanh(x / Xk - 2.0) - TANH_M2) / (1.0 - TANH_M2)


def gp_buildup_factor(x: float, b: float, c: float, a: float, Xk: float, d: float) -> float:
    """B(E, x) for a single element/material whose GP coefficients are known."""
    if x <= 0:
        return 1.0
    K = gp_K(x, c, a, d, Xk)
    if abs(K - 1.0) < 1e-12:
        return 1.0 + (b - 1.0) * x
    return 1.0 + (b - 1.0) * (K**x - 1.0) / (K - 1.0)


def interpolate_gp_coeffs(z1: float, z2: float, zeq: float, coeffs1, coeffs2) -> tuple:
    """Log(Z) interpolation of (b, c, a, Xk, d) between two bracketing elements."""
    if z1 == z2:
        return tuple(coeffs1)
    log_z1 = math.log10(z1)
    log_z2 = math.log10(z2)
    log_zeq = math.log10(zeq)
    frac = (log_zeq - log_z1) / (log_z2 - log_z1)
    return tuple(p1 + (p2 - p1) * frac for p1, p2 in zip(coeffs1, coeffs2))


def element_gp_coeffs(z: int, energy_mev: float, kind: str) -> tuple:
    """(b, c, a, Xk, d) for a single standard element/mixture at an exact tabulated energy.

    `z` is an atomic number (int, for one of the 23 standard elements) or one
    of the mixture keys "WATER" / "AIR" / "CONCRETE". `energy_mev` must be an
    exact key present in GP_COEFFS[kind][z] -- use material_gp_coeffs() below
    for energy interpolation.
    """
    return GP_COEFFS[kind][z][energy_mev]


def _bracket(sorted_values: list[float], target: float) -> tuple[float, float]:
    """Return the two values in a sorted list bracketing target (clamped at the ends)."""
    if target <= sorted_values[0]:
        return sorted_values[0], sorted_values[0]
    if target >= sorted_values[-1]:
        return sorted_values[-1], sorted_values[-1]
    for lo, hi in zip(sorted_values, sorted_values[1:]):
        if lo <= target <= hi:
            return lo, hi
    raise AssertionError("unreachable")


def energy_interpolated_coeffs(z: int, energy_mev: float, kind: str) -> tuple:
    """(b, c, a, Xk, d) for one standard element/mixture at an arbitrary energy,
    log-log interpolated between the two bracketing tabulated energies."""
    table = GP_COEFFS[kind][z]
    energies = sorted(table)
    e_lo, e_hi = _bracket(energies, energy_mev)
    if e_lo == e_hi:
        return table[e_lo]
    c_lo, c_hi = table[e_lo], table[e_hi]
    log_lo, log_hi, log_e = math.log(e_lo), math.log(e_hi), math.log(energy_mev)
    frac = (log_e - log_lo) / (log_hi - log_lo)
    return tuple(v_lo + (v_hi - v_lo) * frac for v_lo, v_hi in zip(c_lo, c_hi))


def material_gp_coeffs(energy_mev: float, zeq: float, z1: int, z2: int, kind: str) -> tuple:
    """Interpolated (b, c, a, Xk, d) for an arbitrary material, given its Zeq(E)
    and the two standard elements (z1 <= zeq <= z2) that bracket it."""
    coeffs1 = energy_interpolated_coeffs(z1, energy_mev, kind)
    coeffs2 = energy_interpolated_coeffs(z2, energy_mev, kind)
    return interpolate_gp_coeffs(z1, z2, zeq, coeffs1, coeffs2)


def buildup_factor(energy_mev: float, x: float, zeq: float, z1: int, z2: int, kind: str) -> float:
    """B(E, x) for an arbitrary material, given its Zeq(E) and GP-coefficient bracket."""
    b, c, a, Xk, d = material_gp_coeffs(energy_mev, zeq, z1, z2, kind)
    return gp_buildup_factor(x, b, c, a, Xk, d)


def element_buildup_factor(z: int, energy_mev: float, x: float, kind: str) -> float:
    """B(E, x) directly for one of the 23 standard elements/3 mixtures (no Zeq
    interpolation needed -- e.g. for a pure-element shield like Lead)."""
    b, c, a, Xk, d = energy_interpolated_coeffs(z, energy_mev, kind)
    return gp_buildup_factor(x, b, c, a, Xk, d)


def buildup_table(zeq_table: dict, mfps: list, kind: str) -> dict:
    """B(E, x) for every (energy, mfp) pair, given a {E: (zeq, z1, z2)} table."""
    results = {}
    for energy, (zeq, z1, z2) in zeq_table.items():
        results[energy] = {
            mfp: buildup_factor(energy, mfp, zeq, z1, z2, kind) for mfp in mfps
        }
    return results
