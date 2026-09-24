"""Python wrapper around NIST's real XCOM v3.1 Fortran program.

XCOM (M.J. Berger & J.H. Hubbell, NIST, 1999) computes photon mass
attenuation coefficients for any element, compound, or mixture at any list
of energies, from its own physics database (vendored in ../vendor/xcom/,
compiled to ../vendor/xcom/xcom -- see vendor/xcom/BUILD.md). This module
drives that binary the same way a human would at the terminal: it feeds the
same stdin prompt sequence XCOM's SPEC subroutine expects, and parses the
cross-section table it writes back out.

Verified against well-known NIST reference points: Pb at 1.0 MeV gives
0.07102 cm^2/g total attenuation (with coherent scattering), water at 1.0
MeV gives 0.07072 cm^2/g -- both exactly the published NIST values.

Protocol notes (reverse-engineered from vendor/xcom/XCOM.f's SPEC
subroutine, since XCOM's own docs don't spell out the stdin sequence):

  - Composition entry has 4 modes (NSUB): 1=element by Z, 2=element by
    symbol, 3=compound by chemical formula, 4=mixture by weight fraction.
    This module only ever drives modes 3 and 4 -- a bare element symbol
    (e.g. "Pb") parses fine as a 1-atom "formula" via mode 3, so there is
    no need for mode 1/2's extra complexity.
  - The "output quantities" prompt (NFORM: barns/atom vs cm2/g) is ONLY
    asked for modes 1/2 in the Fortran source; modes 3/4 skip it and the
    variable defaults to NFORM=3 ("partial interaction coefficients AND
    total attenuation coefficients, all in cm2/g") -- which is exactly the
    mass-attenuation-coefficient quantity this module wants, and for
    modes 3/4 it is applied automatically with no prompt to answer.
  - Energy list entry: NELL=3 ("additional energies only", i.e. skip
    XCOM's built-in standard grid and use exactly the caller's own energy
    list) + INEN=2 ("entry from prepared input file") + a path to a file
    formatted as XCOM's Readme.txt describes: first the count, then the
    energies in MeV, whitespace-separated. XCOM sorts this list itself.
  - Must run with cwd = vendor/xcom/, since the compiled binary opens its
    per-element MDATX3.NNN data files by bare relative name.
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

VENDOR_DIR = Path(__file__).resolve().parents[1] / "vendor" / "xcom"
XCOM_BIN = VENDOR_DIR / "xcom"

# Standard element symbols, Z=1..100 -- only used to let callers pass an
# atomic number instead of a symbol; XCOM itself parses the symbol/formula
# string (not the Z integer) via its own FORM subroutine and internal
# HASH1.DAT/HASH2.DAT symbol tables, so this table's only job is the
# int -> symbol conversion before handing a string to XCOM.
SYMBOL_BY_Z = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
    23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu",
    30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
    44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La",
    58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
    65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
    72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt",
    79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At",
    86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",
    93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es",
    100: "Fm",
}

_NUM = r"([+-]?\d\.\d+E[+-]\d{2})"
_ROW_RE = re.compile(r"^\s*" + _NUM + (r"\s+" + _NUM) * 7 + r"\s*$")

_FIELDS = (
    "energy_MeV",
    "coherent",
    "incoherent",
    "photoelectric",
    "pair_nuclear",
    "pair_electron",
    "total_with_coherent",
    "total_without_coherent",
)


class XcomError(RuntimeError):
    """Raised when the XCOM binary fails or its output can't be parsed."""


def _build_stdin(name: str, composition: str | dict[str, float]) -> list[str]:
    lines = [name]
    if isinstance(composition, dict):
        lines.append("4")  # mixture of elements and/or compounds
        lines.append(str(len(composition)))
        for formula, weight_fraction in composition.items():
            lines.append(formula)
            lines.append(repr(float(weight_fraction)))
        lines.append("1")  # accept, normalize fractions to sum to 1
        # NFORM prompt is skipped for mixtures (defaults to 3).
    else:
        lines.append("3")  # compound, specified by chemical formula
        lines.append(str(composition))
        # NFORM prompt is skipped for formulas (defaults to 3).
    return lines


def _run_xcom(stdin_lines: list[str], energies_mev: list[float]) -> str:
    if not XCOM_BIN.exists():
        raise XcomError(
            f"XCOM binary not found at {XCOM_BIN} -- see vendor/xcom/BUILD.md "
            "to compile it (gfortran -w -std=legacy -o xcom XCOM.f)."
        )
    sorted_energies = sorted(set(float(e) for e in energies_mev))
    with tempfile.NamedTemporaryFile(
        "w", suffix=".txt", delete=False
    ) as energy_f, tempfile.NamedTemporaryFile(
        "w", suffix=".txt", delete=False
    ) as out_f:
        energy_path = Path(energy_f.name)
        out_path = Path(out_f.name)
        energy_f.write(f"{len(sorted_energies)}\n")
        energy_f.write(" ".join(repr(e) for e in sorted_energies) + "\n")

    try:
        stdin_lines = stdin_lines + [
            "3",  # additional energies only (skip XCOM's standard grid)
            "2",  # entry from prepared input file
            str(energy_path),
            str(out_path),
            "1",  # no more output
        ]
        proc = subprocess.run(
            [str(XCOM_BIN)],
            input="\n".join(stdin_lines) + "\n",
            cwd=VENDOR_DIR,
            text=True,
            capture_output=True,
            timeout=30,
        )
        if proc.returncode != 0:
            raise XcomError(
                f"xcom exited with code {proc.returncode}.\n"
                f"stdin was:\n{chr(10).join(stdin_lines)}\n\n"
                f"stdout/stderr:\n{proc.stdout}\n{proc.stderr}"
            )
        if not out_path.exists():
            raise XcomError(
                f"xcom did not produce an output file.\nstdout:\n{proc.stdout}"
            )
        return out_path.read_text()
    finally:
        energy_path.unlink(missing_ok=True)
        out_path.unlink(missing_ok=True)


def _parse_output(text: str) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {field: [] for field in _FIELDS}
    for line in text.splitlines():
        m = _ROW_RE.match(line)
        if not m:
            continue
        values = [float(g) for g in m.groups()]
        for field, value in zip(_FIELDS, values):
            result[field].append(value)
    if not result["energy_MeV"]:
        raise XcomError(f"No data rows parsed from XCOM output:\n{text}")
    return result


def mass_attenuation(
    composition: int | str | dict[str, float], energies_mev: list[float]
) -> dict[str, list[float]]:
    """Mass attenuation coefficients (cm^2/g) for an element, compound, or
    weight-fraction mixture, at the given energies (MeV), via real XCOM.

    composition:
      - int: atomic number Z (converted to its element symbol).
      - str: element symbol ("Pb") or chemical formula ("H2O", "SiO2").
      - dict[str, float]: {formula_or_symbol: weight_fraction, ...} for a
        mixture; fractions need not already sum to 1 (XCOM normalizes them).

    Returns a dict keyed by energy_MeV, coherent, incoherent, photoelectric,
    pair_nuclear, pair_electron, total_with_coherent, total_without_coherent
    -- each a list of floats (cm^2/g, except energy_MeV) aligned by index
    and sorted ascending by energy. total_with_coherent is the standard
    (mu/rho) mass attenuation coefficient; total_without_coherent excludes
    coherent (Rayleigh) scattering.
    """
    if isinstance(composition, int):
        if composition not in SYMBOL_BY_Z:
            raise ValueError(f"Unsupported atomic number: {composition}")
        composition = SYMBOL_BY_Z[composition]

    if isinstance(composition, dict):
        resolved = {}
        for key, frac in composition.items():
            if isinstance(key, int):
                key = SYMBOL_BY_Z[key]
            resolved[key] = frac
        composition = resolved

    stdin_lines = _build_stdin("material", composition)
    raw = _run_xcom(stdin_lines, energies_mev)
    return _parse_output(raw)
