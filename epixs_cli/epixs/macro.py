"""Macro-file batch mode for EpiXS-CLI.

Lets you specify a material once in a plain-text file and re-run the same
calculation without retyping it through the interactive menu -- the same
role Geant4 `.mac` files play for this user's other shielding projects
(compare `~/CHON/Mix_1.txt`'s `/det/addElement C 0.5166` style), so the
syntax deliberately mirrors that: `/command/subcommand value`.

Usage: `./run mymaterial.mac`

Directives (one per line; `#` starts a comment; blank lines ignored):

  /material/name <label>                 -- optional, only used in output
  /material/density <g/cm3>               -- required for HVL/TVL/MFP,
                                              Zeff/Ceff, and neutron removal
  /material/formula <formula-or-symbol>   -- a single element/compound, e.g.
                                              "Pb" or "H2O" -- use this OR
                                              addElement lines, not both
  /material/addElement <Symbol> <weight_fraction>
                                           -- repeatable, builds a mixture

  /energies <list>                        -- MeV, comma/space separated
  /depths <list>                          -- mean free paths, for EBF/EABF

  /output/massAttenuation <file.csv>      -- needs /energies
  /output/hvl <file.csv>                  -- needs /energies, /material/density
  /output/zeff <file.csv>                 -- needs /energies, /material/density
  /output/ebf <file.csv>                  -- needs /energies, /depths
  /output/eabf <file.csv>                 -- needs /energies, /depths
  /output/neutron <file.csv>              -- needs /material/density

Any number of /output/ lines may appear; each is computed and written
independently, in the order given. See ~/epixs_cli/examples/*.mac for
worked examples (the CHON project's ten mixtures).
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


class MacroError(ValueError):
    pass


def _parse_float_list(raw: str) -> list[float]:
    return [float(p) for p in raw.replace(",", " ").split()]


def _write_csv(path: str, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"  wrote {path} ({len(rows)} rows)")


def _compute_mass_attenuation(composition, name, energies, path):
    result = xcom_engine.mass_attenuation(composition, energies)
    header = [
        "Energy (MeV)", "Coherent (cm2/g)", "Incoherent (cm2/g)",
        "Photoelectric (cm2/g)", "Pair-nuclear (cm2/g)", "Pair-electron (cm2/g)",
        "Total w/ coherent (cm2/g)", "Total w/o coherent (cm2/g)",
    ]
    rows = [
        [result["energy_MeV"][i], result["coherent"][i], result["incoherent"][i],
         result["photoelectric"][i], result["pair_nuclear"][i], result["pair_electron"][i],
         result["total_with_coherent"][i], result["total_without_coherent"][i]]
        for i in range(len(result["energy_MeV"]))
    ]
    _write_csv(path, header, rows)


def _compute_hvl(composition, name, energies, density, path):
    result = xcom_engine.mass_attenuation(composition, energies)
    header = ["Energy (MeV)", "mu (1/cm)", "HVL (cm)", "TVL (cm)", "MFP (cm)"]
    rows = []
    for i, e in enumerate(result["energy_MeV"]):
        mu = result["total_with_coherent"][i] * density
        rows.append([e, mu, math.log(2) / mu, math.log(10) / mu, 1.0 / mu])
    _write_csv(path, header, rows)


def _compute_zeff(mass_fracs, density, energies, path):
    ne_static = electron_density(mass_fracs, density)
    header = ["Energy (MeV)", "Zeff", "Zeq", "Bracket Z1", "Bracket Z2",
              "Neff (energy-dep., e/g)", "Ceff (S/m)", "Ne (static, e/cm3)"]
    rows = []
    for e in energies:
        z_eff = zeff(mass_fracs, e, xcom_engine)
        z_eq, z1, z2 = zeq(mass_fracs, e, xcom_engine)
        n_eff = effective_electron_density(mass_fracs, e, xcom_engine)
        c_eff = effective_conductivity(mass_fracs, density, e, xcom_engine)
        rows.append([e, z_eff, z_eq, z1, z2, n_eff, c_eff, ne_static])
    _write_csv(path, header, rows)


def _compute_buildup(mass_fracs, energies, depths, kind, path):
    header = ["Energy (MeV)"] + [f"{x:g} mfp" for x in depths]
    rows = []
    for e in energies:
        z_eq, z1, z2 = zeq(mass_fracs, e, xcom_engine)
        rows.append([e] + [buildup_factor(e, x, z_eq, z1, z2, kind) for x in depths])
    _write_csv(path, header, rows)


def _compute_neutron(mass_fracs, density, name, path):
    sigma_r = removal_cross_section(mass_fracs, density)
    header = ["Material", "Density (g/cm3)", "Sigma_R (1/cm)", "HVL (cm)", "MFP (cm)"]
    _write_csv(path, header, [[name, density, sigma_r, neutron_hvl(sigma_r), neutron_mfp(sigma_r)]])


def run_macro_file(path: str) -> None:
    name = None
    density = None
    formula = None
    elements: dict[str, float] = {}
    energies: list[float] | None = None
    depths: list[float] | None = None
    outputs: list[tuple[str, str]] = []

    with open(path) as f:
        for lineno, raw_line in enumerate(f, 1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                raise MacroError(f"{path}:{lineno}: expected '/command value', got {raw_line!r}")
            cmd, value = parts
            try:
                if cmd == "/material/name":
                    name = value
                elif cmd == "/material/density":
                    density = float(value)
                elif cmd == "/material/formula":
                    formula = value
                elif cmd == "/material/addElement":
                    sym, frac = value.split()
                    elements[sym] = float(frac)
                elif cmd == "/energies":
                    energies = _parse_float_list(value)
                elif cmd == "/depths":
                    depths = _parse_float_list(value)
                elif cmd.startswith("/output/"):
                    outputs.append((cmd[len("/output/"):], value))
                else:
                    raise MacroError(f"{path}:{lineno}: unknown directive {cmd!r}")
            except ValueError as exc:
                raise MacroError(f"{path}:{lineno}: {exc}") from exc

    if formula and elements:
        raise MacroError("use either /material/formula or /material/addElement, not both")
    if not formula and not elements:
        raise MacroError("no material specified (need /material/formula or /material/addElement)")

    if formula:
        composition = formula
        mass_fracs = formula_mass_fractions(formula)
        display_name = name or formula
    else:
        composition = elements
        try:
            mass_fracs = expand_mixture(elements)
        except FormulaError as exc:
            raise MacroError(str(exc)) from exc
        display_name = name or "+".join(f"{s}({w:g})" for s, w in elements.items())

    print(f"Material: {display_name}")
    for quantity, out_path in outputs:
        print(f"Computing {quantity} -> {out_path}")
        if quantity == "massAttenuation":
            if energies is None:
                raise MacroError("/output/massAttenuation needs /energies")
            _compute_mass_attenuation(composition, display_name, energies, out_path)
        elif quantity == "hvl":
            if energies is None or density is None:
                raise MacroError("/output/hvl needs /energies and /material/density")
            _compute_hvl(composition, display_name, energies, density, out_path)
        elif quantity == "zeff":
            if energies is None or density is None:
                raise MacroError("/output/zeff needs /energies and /material/density")
            _compute_zeff(mass_fracs, density, energies, out_path)
        elif quantity == "ebf":
            if energies is None or depths is None:
                raise MacroError("/output/ebf needs /energies and /depths")
            _compute_buildup(mass_fracs, energies, depths, "exposure", out_path)
        elif quantity == "eabf":
            if energies is None or depths is None:
                raise MacroError("/output/eabf needs /energies and /depths")
            _compute_buildup(mass_fracs, energies, depths, "absorption", out_path)
        elif quantity == "neutron":
            if density is None:
                raise MacroError("/output/neutron needs /material/density")
            _compute_neutron(mass_fracs, density, display_name, out_path)
        else:
            raise MacroError(f"unknown /output/{quantity}")


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python3 macro.py <file.mac>", file=sys.stderr)
        sys.exit(1)
    try:
        run_macro_file(sys.argv[1])
    except (MacroError, xcom_engine.XcomError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
