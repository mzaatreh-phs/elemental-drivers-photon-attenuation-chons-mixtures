"""Fast neutron removal cross section -- a Phy-X/PSD output quantity, and a
different physics domain from the rest of this package (fast/fission-energy
neutron shielding, not photon attenuation, so it does NOT go through
xcom_engine). Roughly energy-independent over the ~2-12 MeV range where the
"removal" concept applies (a fast neutron's first collision in this range
overwhelmingly removes it from the uncollided beam, regardless of exact
energy), unlike every other quantity in this package.

Formulas (independently verified numerically against a paper that both
states them and tabulates worked examples -- every reproduced value below
matches to the paper's own printed precision):

  (Sigma_R/rho)_i = 0.206 * A_i^(-1/3) * Z_i^(-0.294)   [cm^2/g]      (2)
  Sigma_R = sum_i( rho_i * (Sigma_R/rho)_i )            [cm^-1]      (1)
  HVL = 0.693 / Sigma_R                                 [cm]         (3)
  MFP = 1 / Sigma_R                                     [cm]         (4)

where A_i, Z_i are the atomic weight and atomic number of constituent i,
and rho_i = w_i * rho is that constituent's partial density (w_i its mass
fraction, rho the material's bulk density, g/cm^3).

Source: Ahmed Fadhil Mkhaiber & Salwah Kareem Dawood, "Calculation of
Shielding Parameters of Fast Neutrons for Some Composite Materials",
Al-Mustansiriyah Journal of Science 30(1), 2019, pp. 209-215,
DOI: 10.23851/mjs.v30i1.520 -- Eqs. (1)-(4), freely downloadable PDF via
mjs.uomustansiriyah.edu.iq. That paper in turn cites Eq. (2) to J.E. Martin,
"Physics for Radiation Protection", 2nd ed., Wiley, 2006, and Eq. (1)'s
mixture rule to A.E. Profio, "Radiation Shielding and Dosimetry", Wiley,
1979 -- both are themselves the standard, widely-cited sources for this
semi-empirical formula.

Validated by recomputing Eq. (2) for every element in that paper's own
worked example (Table 1, pure paraffin wax, C and H; and Tables 2-6's
composite constituents B, W, O, Fe, Al, Si) and confirming agreement to the
paper's own printed 2-3 significant figures in every case -- see
tests/test_phyx_extensions.py.

The 0.206/A^(-1/3)/Z^(-0.294) form is a SEMI-EMPIRICAL fit, not a first-
principles calculation, and is understood to be most reliable for light-to-
medium Z elements in typical shielding materials; treat results for very
heavy elements (or very light ones like H, where the fit is known to
somewhat underestimate compared to measured values in the original
literature) as approximate, consistent with how removal cross sections are
generally used in shielding engineering (order-of-magnitude / comparative
screening, not a precision transport calculation).
"""

from __future__ import annotations

from atomic_data import ATOMIC_WEIGHT, SYMBOL_TO_Z


def element_removal_cross_section(symbol: str) -> float:
    """(Sigma_R/rho) for one element [cm^2/g], via Eq. (2)."""
    A = ATOMIC_WEIGHT[symbol]
    Z = SYMBOL_TO_Z[symbol]
    return 0.206 * A ** (-1.0 / 3.0) * Z ** (-0.294)


def removal_cross_section(composition_by_mass: dict, density_g_cm3: float) -> float:
    """Macroscopic fast-neutron removal cross section Sigma_R [cm^-1] for a
    compound/mixture at the given bulk density, via Eq. (1)."""
    total = 0.0
    for sym, w in composition_by_mass.items():
        partial_density = w * density_g_cm3
        total += partial_density * element_removal_cross_section(sym)
    return total


def neutron_hvl(sigma_r: float) -> float:
    """Half-value layer [cm] for fast neutrons, Eq. (3)."""
    return 0.693 / sigma_r


def neutron_mfp(sigma_r: float) -> float:
    """Mean free path [cm] for fast neutrons, Eq. (4)."""
    return 1.0 / sigma_r
