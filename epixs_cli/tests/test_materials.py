import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "epixs"))

import pytest

from materials import FormulaError, expand_mixture, formula_mass_fractions, parse_formula


def test_parse_simple_formula():
    assert parse_formula("H2O") == {"H": 2.0, "O": 1.0}


def test_parse_bare_element():
    assert parse_formula("Pb") == {"Pb": 1.0}


def test_parse_multi_digit_subscript():
    assert parse_formula("Nb2O5") == {"Nb": 2.0, "O": 5.0}


def test_parse_parenthesized_group():
    assert parse_formula("Ca(OH)2") == {"Ca": 1.0, "O": 2.0, "H": 2.0}


def test_parse_unknown_element_raises():
    with pytest.raises(FormulaError):
        parse_formula("Xx2O")


def test_parse_unbalanced_parens_raises():
    with pytest.raises(FormulaError):
        parse_formula("Ca(OH2")


def test_formula_mass_fractions_water_matches_known_value():
    fracs = formula_mass_fractions("H2O")
    assert fracs["H"] == pytest.approx(0.111898, rel=2e-3)
    assert fracs["O"] == pytest.approx(0.888102, rel=2e-3)
    assert sum(fracs.values()) == pytest.approx(1.0)


def test_expand_mixture_normalizes_and_combines_shared_elements():
    combined = expand_mixture({"H2O": 0.9, "NaCl": 0.1})
    assert sum(combined.values()) == pytest.approx(1.0)
    assert set(combined) == {"H", "O", "Na", "Cl"}


def test_expand_mixture_unnormalized_fractions():
    a = expand_mixture({"O": 1, "P": 1})
    b = expand_mixture({"O": 10, "P": 10})
    assert a == pytest.approx(b)
