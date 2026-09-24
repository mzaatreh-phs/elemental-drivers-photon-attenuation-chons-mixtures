"""
Equivalent atomic number (Zeq), effective atomic number (Zeff), and electron
density for an arbitrary compound/mixture, at a given photon energy.

Formulas verified against an open-access, citable source (not taken from
memory):

  Kaewkhao/... et al., "Radiation Shielding Properties and Exposure Buildup
  Factor of TI-Al-Nb Alloy Materials", IOSR Journal of Applied Physics,
  Vol. 11, Issue 2, Series I (2019), pp. 68-74, DOI: 10.9790/4861-1102016874,
  https://www.iosrjournals.org/iosr-jap/papers/Vol11-issue2/Series-1/
  K1102016874.pdf -- Eqs. (1)-(8) transcribed 2026-09-10 (pdftotext of the
  open-access PDF; this exact set of equations, with the same symbols, is
  reproduced near-verbatim across dozens of papers in this field, e.g. the
  Sidhu et al. Zeq formula and the Manohara & Hanagodimath Zeff formula).

  (1) (mu/rho)_mix = sum_i w_i (mu/rho)_i                    [Bragg additivity]
  (2) sigma_t,a  = (mu/rho)_mix / (NA * sum_i w_i/A_i)        [total atomic cross section, cm^2/atom]
  (3) sigma_t,el = (1/NA) * sum_i (f_i A_i / Z_i) (mu/rho)_i  [total electronic cross section, cm^2/electron]
  (4) Zeff = sigma_t,a / sigma_t,el
  (5) Nel  = (mu/rho)_mix / sigma_t,el                        [effective electron density, electrons/g]
  (7) Zeq  = [Z1(log R2 - log R) + Z2(log R - log R1)] / (log R2 - log R1)
       where R = (mu/rho)_Compton / (mu/rho)_Total, R1/R2 the same ratio for
       the two standard elements Z1 <= Zeq <= Z2 bracketing R.

w_i is weight (mass) fraction, f_i is mole (atom) fraction, A_i atomic
weight, Z_i atomic number, NA Avogadro's number -- all per constituent
element i. Zeff and Nel are ENERGY-DEPENDENT (they vary with photon energy
because the interaction cross sections do); this is a different quantity
from the plain physical electron density (Ne, below), which is a static
mass-composition fact independent of photon energy.

Zeff and Nel formulas were independently re-derived from Eqs. (1)-(5) above
to confirm self-consistency (numerator/denominator both reduce cleanly to
NA*sigma_t,a and NA*sigma_t,el respectively) -- see the derivation notes
kept alongside this project's session records if you need to re-verify.

This module needs, for each constituent element, the per-element mass
attenuation coefficients (mu/rho)_i, split into "incoherent" (Compton) and
"total", at the energy of interest. That is the job of the XCOM engine
(epixs/xcom_engine.py) being built in parallel elsewhere in this project;
if it is not importable yet, the functions here that need it raise
ImportError with a clear message rather than failing silently or
fabricating numbers.
"""

from __future__ import annotations

import math

from atomic_data import ATOMIC_WEIGHT
from gp_data import GP_COEFFS, STANDARD_ELEMENTS_Z

AVOGADRO = 6.02214076e23  # mol^-1


def _get_xcom_engine():
    try:
        import xcom_engine  # noqa: F401  (built in parallel; may not exist yet)
    except ImportError as exc:
        raise ImportError(
            "zeff.py needs epixs.xcom_engine.mass_attenuation() for per-element "
            "mu/rho lookups (Zeq/Zeff/Nel all depend on it) -- that module is "
            "being built separately and isn't importable yet."
        ) from exc
    return xcom_engine


def mole_fractions(composition_by_mass: dict) -> dict:
    """{symbol: weight_fraction} -> {symbol: mole_fraction}."""
    moles = {sym: w / ATOMIC_WEIGHT[sym] for sym, w in composition_by_mass.items()}
    total = sum(moles.values())
    return {sym: n / total for sym, n in moles.items()}


def mean_atomic_weight(composition_by_mass: dict) -> float:
    """<A> = 1 / sum_i(w_i / A_i) -- the mean molar mass per atom of the mixture."""
    return 1.0 / sum(w / ATOMIC_WEIGHT[sym] for sym, w in composition_by_mass.items())


def electron_density(composition_by_mass: dict, density_g_cm3: float) -> float:
    """Plain physical electron density Ne [electrons/cm^3] -- energy-independent,
    a pure mass-composition fact: Ne = rho * NA * sum_i(w_i * Z_i / A_i)."""
    from atomic_data import SYMBOL_TO_Z

    per_gram = AVOGADRO * sum(
        w * SYMBOL_TO_Z[sym] / ATOMIC_WEIGHT[sym] for sym, w in composition_by_mass.items()
    )
    return per_gram * density_g_cm3


def _mixture_mu_rho(composition_by_mass: dict, energy_mev: float, component: str, xcom_engine) -> float:
    """(mu/rho)_mix for one XCOM component ('incoherent' or 'total') via Bragg additivity, Eq. (1)."""
    total = 0.0
    for sym, w in composition_by_mass.items():
        per_element = xcom_engine.mass_attenuation(sym, [energy_mev])
        total += w * per_element[component][0]
    return total


def compton_to_total_ratio(composition_by_mass: dict, energy_mev: float, xcom_engine=None) -> float:
    """R = (mu/rho)_Compton / (mu/rho)_Total for a compound/mixture at one energy.

    "Total" here means total attenuation WITHOUT coherent (Rayleigh)
    scattering -- coherent scattering is elastic (no energy transfer, no
    change in photon population/energy), so the standard Zeq convention
    excludes it from this ratio's denominator. Confirmed empirically against
    a real EpiXS export (CHON project's MIX01, 30 keV): using
    total_with_coherent gives R=0.5784 vs EpiXS's own printed R=0.66055;
    using total_without_coherent gives R=0.66060 -- matching to 4 sig figs.
    (An earlier version of this function used total_with_coherent, which
    reproduced EpiXS's Lead/water single-element buildup factors fine only
    because coherent scattering happens to be negligible there at the
    energies tested, but caused up to ~35% buildup-factor error for real
    low-Z mixtures at low energy, where coherent scattering is not
    negligible relative to incoherent.)
    """
    xcom_engine = xcom_engine or _get_xcom_engine()
    incoh = _mixture_mu_rho(composition_by_mass, energy_mev, "incoherent", xcom_engine)
    total = _mixture_mu_rho(composition_by_mass, energy_mev, "total_without_coherent", xcom_engine)
    return incoh / total


def element_compton_to_total_ratio(z: int, energy_mev: float, xcom_engine=None) -> float:
    """R for a single standard element (used to build the R1/R2 bracket for Zeq).
    See compton_to_total_ratio() docstring for why the denominator excludes
    coherent scattering."""
    xcom_engine = xcom_engine or _get_xcom_engine()
    from atomic_data import ELEMENTS

    sym = ELEMENTS[z][0]
    per_element = xcom_engine.mass_attenuation(sym, [energy_mev])
    return per_element["incoherent"][0] / per_element["total_without_coherent"][0]


def zeq(composition_by_mass: dict, energy_mev: float, xcom_engine=None) -> tuple[float, int, int]:
    """Equivalent atomic number Zeq(E), Eq. (7), plus the (Z1, Z2) bracket it
    was interpolated within (Z1 <= Zeq <= Z2, both in STANDARD_ELEMENTS_Z)."""
    xcom_engine = xcom_engine or _get_xcom_engine()
    R = compton_to_total_ratio(composition_by_mass, energy_mev, xcom_engine)

    ratios = [
        (z, element_compton_to_total_ratio(z, energy_mev, xcom_engine))
        for z in STANDARD_ELEMENTS_Z
    ]
    # R decreases monotonically with Z at fixed energy (Compton fraction falls
    # as photoelectric/pair-production grow with Z) -- sort ascending in Z, so
    # r is descending along this list; the bracket search below is written to
    # match that direction (NOT sorted-by-R, which would silently pair each
    # z1 with a z2 < z1 and feed the interpolation formula an inverted bracket).
    ratios.sort(key=lambda zr: zr[0])

    if R >= ratios[0][1]:
        z_lo, r_lo = ratios[0]
        return float(z_lo), z_lo, z_lo
    if R <= ratios[-1][1]:
        z_hi, r_hi = ratios[-1]
        return float(z_hi), z_hi, z_hi

    for (z1, r1), (z2, r2) in zip(ratios, ratios[1:]):
        if r2 <= R <= r1:
            if r1 == r2:
                return float(z1), z1, z2
            log_r1, log_r2, log_r = math.log10(r1), math.log10(r2), math.log10(R)
            z_eq = (z1 * (log_r2 - log_r) + z2 * (log_r - log_r1)) / (log_r2 - log_r1)
            return z_eq, z1, z2

    raise AssertionError("unreachable: R not bracketed despite range checks above")


def zeff(composition_by_mass: dict, energy_mev: float, xcom_engine=None) -> float:
    """Effective atomic number Zeff(E), Eq. (4): sigma_t,a / sigma_t,el."""
    xcom_engine = xcom_engine or _get_xcom_engine()
    from atomic_data import SYMBOL_TO_Z

    f = mole_fractions(composition_by_mass)
    numerator = 0.0  # ~ sigma_t,a * NA
    denominator = 0.0  # ~ sigma_t,el * NA
    for sym, w in composition_by_mass.items():
        mu_rho_i = xcom_engine.mass_attenuation(sym, [energy_mev])["total_with_coherent"][0]
        A_i = ATOMIC_WEIGHT[sym]
        Z_i = SYMBOL_TO_Z[sym]
        numerator += f[sym] * A_i * mu_rho_i
        denominator += f[sym] * (A_i / Z_i) * mu_rho_i
    return numerator / denominator


def effective_electron_density(composition_by_mass: dict, energy_mev: float, xcom_engine=None) -> float:
    """Nel(E), Eq. (5): effective electron density [electrons/g] -- energy
    dependent, distinct from the plain physical electron_density() above."""
    xcom_engine = xcom_engine or _get_xcom_engine()
    from atomic_data import SYMBOL_TO_Z

    f = mole_fractions(composition_by_mass)
    mu_rho_mix = _mixture_mu_rho(composition_by_mass, energy_mev, "total_with_coherent", xcom_engine)
    sigma_t_el_times_NA = 0.0
    for sym, w in composition_by_mass.items():
        mu_rho_i = xcom_engine.mass_attenuation(sym, [energy_mev])["total_with_coherent"][0]
        A_i = ATOMIC_WEIGHT[sym]
        Z_i = SYMBOL_TO_Z[sym]
        sigma_t_el_times_NA += f[sym] * (A_i / Z_i) * mu_rho_i
    return mu_rho_mix / sigma_t_el_times_NA * AVOGADRO


# Physical constants for effective_conductivity() (SI).
_ELECTRON_CHARGE_C = 1.602176634e-19
_ELECTRON_MASS_KG = 9.1093837015e-31
_PLANCK_H_JS = 6.62607015e-34
_BOLTZMANN_K_JPERK = 1.380649e-23
_ROOM_TEMPERATURE_K = 300.0


def effective_conductivity(composition_by_mass: dict, density_g_cm3: float,
                            energy_mev: float, xcom_engine=None) -> float:
    """Ceff(E) -- "effective conductivity", a Phy-X/PSD output quantity.

    Reproduces, exactly as published, the pair of equations used across the
    Phy-X/PSD-adjacent shielding literature (independently cross-checked
    against two sources that state the same formula):

      tau = h / (2*pi * k * T), T = 300 K            [electron relaxation
                                                        time "on the Fermi
                                                        surface", evaluated
                                                        at room temperature]
      Ceff = (Neff * rho * e^2 * tau / m_e) * 1e3

    where Neff is THIS module's effective_electron_density() (electrons/g,
    energy-dependent), rho is the material density (g/cm^3), e and m_e are
    the electron charge and rest mass, and h/k are Planck's/Boltzmann's
    constants.

    Sources (independently found and cross-checked -- both state the same
    Ceff equation in the same symbols; the tau formula came from the second,
    which gives it explicitly as tau = h/(600*pi*k), i.e. h/(2*pi*k*300)):
      - Al-Buriahi et al. (2025), "Assessment of gamma and neutron shielding
        features of Se80-yTe20My alloys", Sci. Rep. 15, 34075, Eq. (9a)-(9b)
        (nature.com/articles/s41598-025-34075-3).
      - Rammah, Y.S. et al. (2021), "Trivalent Ions and Their Impacts on
        Effective Conductivity at 300 K ... Bismo-Borate Glasses", Materials
        14(19), 5894, Eq. (6)-(7) (PMC8510023).

    NOT independently re-derived or unit-checked against SI conductivity
    (S/m) from first principles here -- this function reproduces the
    published formula's numbers as-is, per this project's policy of citing
    rather than re-deriving specialized literature quantities. No local
    ground-truth export (e.g. from EpiXS) was available to validate this
    one against, unlike Zeff/Neff/Zeq above -- treat with proportionally
    more caution than those.
    """
    n_eff = effective_electron_density(composition_by_mass, energy_mev, xcom_engine)
    tau = _PLANCK_H_JS / (2 * math.pi * _BOLTZMANN_K_JPERK * _ROOM_TEMPERATURE_K)
    return (n_eff * density_g_cm3 * _ELECTRON_CHARGE_C**2 * tau / _ELECTRON_MASS_KG) * 1e3
