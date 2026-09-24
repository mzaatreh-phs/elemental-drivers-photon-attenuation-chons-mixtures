#!/usr/bin/env python3
"""EpiXS-CLI: a terminal-native gamma-ray shielding calculator.

A guided-menu replacement for EpiXS/EPICS2017 (Windows GUI tools this
project's user previously ran under Wine): mass attenuation coefficients,
effective atomic number / electron density, and G-P exposure/energy-
absorption buildup factors, for elements, compounds, and weight-fraction
mixtures, computed from the real NIST XCOM engine (xcom_engine.py) and the
digitized ANSI/ANS-6.4.3-1991 G-P buildup table (gp_data.py) -- see those
modules' docstrings for sourcing and validation details.
"""

from __future__ import annotations

import csv
import math
import sys

import xcom_engine
from gp_calculator import buildup_factor
from materials import FormulaError, expand_mixture, formula_mass_fractions
from neutron import neutron_hvl, neutron_mfp, removal_cross_section
from zeff import effective_conductivity, effective_electron_density, electron_density, zeff, zeq

BANNER = r"""
==============================================================
  EpiXS-CLI -- gamma-ray attenuation & buildup factor toolkit
==============================================================
"""


def _prompt(msg: str) -> str:
    try:
        return input(msg).strip()
    except EOFError:
        print()
        sys.exit(0)


def _prompt_float_list(msg: str) -> list[float]:
    raw = _prompt(msg)
    parts = raw.replace(",", " ").split()
    try:
        return [float(p) for p in parts]
    except ValueError as exc:
        print(f"  Could not parse a number in {raw!r}: {exc}")
        return _prompt_float_list(msg)


def _prompt_choice(msg: str, choices: list[str]) -> int:
    print(msg)
    for i, choice in enumerate(choices, 1):
        print(f"    {i}. {choice}")
    raw = _prompt("  Enter choice: ")
    try:
        n = int(raw)
        if 1 <= n <= len(choices):
            return n
    except ValueError:
        pass
    print("  Invalid choice, try again.")
    return _prompt_choice(msg, choices)


def _prompt_material() -> tuple[str | dict, dict[str, float], str]:
    """Returns (xcom_composition, composition_by_mass, display_name)."""
    kind = _prompt_choice(
        "\nHow is the material specified?",
        ["Element (e.g. Pb)", "Compound, chemical formula (e.g. H2O)", "Mixture, by weight fraction"],
    )
    if kind in (1, 2):
        formula = _prompt("  Enter symbol/formula: ")
        try:
            mass_fracs = formula_mass_fractions(formula)
        except FormulaError as exc:
            print(f"  {exc}")
            return _prompt_material()
        return formula, mass_fracs, formula

    n = _prompt("  How many components? ")
    try:
        n = int(n)
    except ValueError:
        print("  Enter an integer.")
        return _prompt_material()

    components: dict[str, float] = {}
    for i in range(1, n + 1):
        formula = _prompt(f"  Component {i} symbol/formula: ")
        frac = _prompt(f"  Component {i} weight fraction: ")
        try:
            components[formula] = float(frac)
        except ValueError:
            print("  Enter a number.")
            return _prompt_material()
    try:
        mass_fracs = expand_mixture(components)
    except FormulaError as exc:
        print(f"  {exc}")
        return _prompt_material()
    name = "+".join(f"{f}({w:g})" for f, w in components.items())
    return components, mass_fracs, name


def _maybe_save_csv(header: list[str], rows: list[list]) -> None:
    save = _prompt("\nSave results to CSV? (y/N): ").lower()
    if save != "y":
        return
    path = _prompt("  File name: ")
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"  Wrote {path}")


def run_mass_attenuation() -> None:
    composition, _mass_fracs, name = _prompt_material()
    energies = _prompt_float_list("\nEnergies (MeV), space/comma separated: ")
    result = xcom_engine.mass_attenuation(composition, energies)

    header = [
        "Energy (MeV)", "Coherent (cm2/g)", "Incoherent (cm2/g)",
        "Photoelectric (cm2/g)", "Pair-nuclear (cm2/g)", "Pair-electron (cm2/g)",
        "Total w/ coherent (cm2/g)", "Total w/o coherent (cm2/g)",
    ]
    print(f"\n=== Mass attenuation coefficients: {name} ===")
    print(f"{'E (MeV)':>10} {'Coherent':>10} {'Incoher.':>10} {'Photoel.':>10} "
          f"{'PairNuc':>10} {'PairEl':>10} {'Tot(w/coh)':>12} {'Tot(no coh)':>12}")
    rows = []
    for i, e in enumerate(result["energy_MeV"]):
        row = [
            e, result["coherent"][i], result["incoherent"][i], result["photoelectric"][i],
            result["pair_nuclear"][i], result["pair_electron"][i],
            result["total_with_coherent"][i], result["total_without_coherent"][i],
        ]
        rows.append(row)
        print(f"{e:10.4g} {row[1]:10.4g} {row[2]:10.4g} {row[3]:10.4g} "
              f"{row[4]:10.4g} {row[5]:10.4g} {row[6]:12.5g} {row[7]:12.5g}")
    _maybe_save_csv(header, rows)

    want_hvl = _prompt("\nAlso compute linear attenuation coeff / HVL / TVL / MFP? (y/N): ").lower()
    if want_hvl == "y":
        density = None
        while density is None:
            raw = _prompt("  Material density (g/cm^3): ")
            try:
                density = float(raw)
            except ValueError:
                print("  Enter a number.")
        hvl_header = ["Energy (MeV)", "mu (1/cm)", "HVL (cm)", "TVL (cm)", "MFP (cm)"]
        print(f"\n{'E (MeV)':>10} {'mu (1/cm)':>12} {'HVL (cm)':>12} {'TVL (cm)':>12} {'MFP (cm)':>12}")
        hvl_rows = []
        for i, e in enumerate(result["energy_MeV"]):
            mu = result["total_with_coherent"][i] * density
            hvl, tvl, mfp = math.log(2) / mu, math.log(10) / mu, 1.0 / mu
            hvl_rows.append([e, mu, hvl, tvl, mfp])
            print(f"{e:10.4g} {mu:12.5g} {hvl:12.5g} {tvl:12.5g} {mfp:12.5g}")
        _maybe_save_csv(hvl_header, hvl_rows)


def run_zeff_ne() -> None:
    _composition, mass_fracs, name = _prompt_material()
    density = None
    while density is None:
        raw = _prompt("\nMaterial density (g/cm^3), for physical electron density: ")
        try:
            density = float(raw)
        except ValueError:
            print("  Enter a number.")
    energies = _prompt_float_list("Energies (MeV), space/comma separated: ")

    ne_static = electron_density(mass_fracs, density)
    header = ["Energy (MeV)", "Zeff", "Zeq", "Bracket Z1", "Bracket Z2",
               "Neff (energy-dep., e/g)", "Ceff (S/m)", "Ne (static, e/cm3)"]
    print(f"\n=== Effective atomic number / electron density: {name} ===")
    print(f"  Static electron density Ne = {ne_static:.6g} electrons/cm^3 (energy-independent)")
    print(f"{'E (MeV)':>10} {'Zeff':>8} {'Zeq':>8} {'Z1':>5} {'Z2':>5} {'Neff (e/g)':>14} {'Ceff (S/m)':>14}")
    rows = []
    for e in energies:
        z_eff = zeff(mass_fracs, e, xcom_engine)
        z_eq, z1, z2 = zeq(mass_fracs, e, xcom_engine)
        n_eff = effective_electron_density(mass_fracs, e, xcom_engine)
        c_eff = effective_conductivity(mass_fracs, density, e, xcom_engine)
        rows.append([e, z_eff, z_eq, z1, z2, n_eff, c_eff, ne_static])
        print(f"{e:10.4g} {z_eff:8.3f} {z_eq:8.3f} {z1:5d} {z2:5d} {n_eff:14.5g} {c_eff:14.5g}")
    _maybe_save_csv(header, rows)


def run_buildup_factor(kind: str, label: str) -> None:
    _composition, mass_fracs, name = _prompt_material()
    energies = _prompt_float_list("\nEnergies (MeV), space/comma separated: ")
    mfps = _prompt_float_list("Penetration depths (mean free paths), space/comma separated: ")

    header = ["Energy (MeV)"] + [f"{x:g} mfp" for x in mfps]
    print(f"\n=== {label}: {name} ===")
    col_w = 12
    print(f"{'E (MeV)':>10}" + "".join(f"{x:g} mfp".rjust(col_w) for x in mfps))
    rows = []
    for e in energies:
        z_eq, z1, z2 = zeq(mass_fracs, e, xcom_engine)
        row_vals = [buildup_factor(e, x, z_eq, z1, z2, kind) for x in mfps]
        rows.append([e] + row_vals)
        print(f"{e:10.4g}" + "".join(f"{v:.4g}".rjust(col_w) for v in row_vals))
    _maybe_save_csv(header, rows)


def run_neutron_removal() -> None:
    _composition, mass_fracs, name = _prompt_material()
    density = None
    while density is None:
        raw = _prompt("\nMaterial density (g/cm^3): ")
        try:
            density = float(raw)
        except ValueError:
            print("  Enter a number.")

    sigma_r = removal_cross_section(mass_fracs, density)
    hvl = neutron_hvl(sigma_r)
    mfp = neutron_mfp(sigma_r)
    print(f"\n=== Fast neutron removal cross section: {name} ===")
    print("  (roughly energy-independent, valid for the ~2-12 MeV fast-neutron range)")
    print(f"  Sigma_R (macroscopic removal cross section) = {sigma_r:.5g} cm^-1")
    print(f"  Half-value layer (HVL)                       = {hvl:.5g} cm")
    print(f"  Mean free path (MFP)                          = {mfp:.5g} cm")
    _maybe_save_csv(
        ["Material", "Density (g/cm3)", "Sigma_R (1/cm)", "HVL (cm)", "MFP (cm)"],
        [[name, density, sigma_r, hvl, mfp]],
    )


def main() -> None:
    print(BANNER)
    while True:
        choice = _prompt_choice(
            "Main menu -- what would you like to compute?",
            [
                "Mass attenuation coefficient",
                "Effective atomic number (Zeff), electron density & conductivity",
                "Exposure buildup factor (EBF)",
                "Energy absorption buildup factor (EABF)",
                "Fast neutron removal cross section",
                "Exit",
            ],
        )
        try:
            if choice == 1:
                run_mass_attenuation()
            elif choice == 2:
                run_zeff_ne()
            elif choice == 3:
                run_buildup_factor("exposure", "Exposure Buildup Factor (EBF)")
            elif choice == 4:
                run_buildup_factor("absorption", "Energy Absorption Buildup Factor (EABF)")
            elif choice == 5:
                run_neutron_removal()
            elif choice == 6:
                print("Goodbye.")
                return
        except xcom_engine.XcomError as exc:
            print(f"\n  XCOM error: {exc}")


if __name__ == "__main__":
    main()
