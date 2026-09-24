"""Chemical formula parsing and mixture composition helpers.

XCOM's own Fortran FORM subroutine parses formula strings internally (see
xcom_engine.py), so mass_attenuation() never needs this module. But zeff.py's
Zeff/Zeq/electron-density functions need an explicit {element_symbol:
weight_fraction} dict, so the CLI needs its own formula -> elemental
weight-fraction expansion to build that dict from user input like "H2O" or
a mixture of formulas.
"""

from __future__ import annotations

import re

from atomic_data import ATOMIC_WEIGHT, SYMBOL_TO_Z

_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)|(\()|(\))(\d*\.?\d*)")


class FormulaError(ValueError):
    pass


def parse_formula(formula: str) -> dict[str, float]:
    """Chemical formula (e.g. "H2O", "SiO2", "Ca(OH)2") -> {symbol: atom_count}.

    Supports simple nesting via one level of parentheses with a trailing
    multiplier, e.g. "Ca(OH)2". Counts default to 1 when omitted.
    """
    formula = formula.strip()
    if not formula:
        raise FormulaError("empty formula")

    stack: list[dict[str, float]] = [{}]
    pos = 0
    while pos < len(formula):
        m = _TOKEN_RE.match(formula, pos)
        if not m or m.start() != pos:
            raise FormulaError(f"could not parse {formula!r} at position {pos}")
        sym, count, open_paren, close_paren, close_count = m.groups()
        if open_paren:
            stack.append({})
        elif close_paren:
            group = stack.pop()
            mult = float(close_count) if close_count else 1.0
            target = stack[-1]
            for s, n in group.items():
                target[s] = target.get(s, 0.0) + n * mult
        elif sym:
            if sym not in SYMBOL_TO_Z:
                raise FormulaError(f"unknown element symbol {sym!r} in {formula!r}")
            n = float(count) if count else 1.0
            target = stack[-1]
            target[sym] = target.get(sym, 0.0) + n
        else:
            raise FormulaError(f"could not parse {formula!r} at position {pos}")
        pos = m.end()

    if len(stack) != 1:
        raise FormulaError(f"unbalanced parentheses in {formula!r}")
    atom_counts = stack[0]
    if not atom_counts:
        raise FormulaError(f"no elements found in {formula!r}")
    return atom_counts


def formula_mass_fractions(formula: str) -> dict[str, float]:
    """Chemical formula -> {element_symbol: mass_fraction}, normalized to sum 1.

    A bare element symbol (e.g. "Pb") is a valid one-atom "formula" here.
    """
    atom_counts = parse_formula(formula)
    masses = {sym: n * ATOMIC_WEIGHT[sym] for sym, n in atom_counts.items()}
    total = sum(masses.values())
    return {sym: m / total for sym, m in masses.items()}


def expand_mixture(components: dict[str, float]) -> dict[str, float]:
    """{formula_or_symbol: weight_fraction} -> flat {element_symbol: weight_fraction},
    normalized to sum 1. Each component's own weight fraction need not already
    sum to 1 across components (normalized here), matching XCOM's own
    "accept, normalize" mixture behavior.
    """
    total_weight = sum(components.values())
    if total_weight <= 0:
        raise FormulaError("mixture weight fractions must sum to a positive value")

    combined: dict[str, float] = {}
    for formula, weight in components.items():
        norm_weight = weight / total_weight
        for sym, frac in formula_mass_fractions(formula).items():
            combined[sym] = combined.get(sym, 0.0) + norm_weight * frac
    return combined
