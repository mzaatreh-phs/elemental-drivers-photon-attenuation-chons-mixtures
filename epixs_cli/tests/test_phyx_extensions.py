"""Tests for the Phy-X/PSD-parity extensions: effective conductivity (Ceff)
and fast neutron removal cross section, added to match Phy-X/PSD's full
parameter list. See epixs/zeff.py's effective_conductivity() and
epixs/neutron.py docstrings for sourcing.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "epixs"))

import pytest

from neutron import element_removal_cross_section, neutron_hvl, neutron_mfp, removal_cross_section

# ---------------------------------------------------------------------------
# Fast neutron removal cross section: validated against Mkhaiber & Dawood
# (2019), Al-Mustansiriyah J. Science 30(1), Table 1 (pure paraffin wax,
# C/H) and the per-element (Sigma_R/rho) values listed alongside Tables 2-6
# (B, W, O, Fe, Al, Si) -- every value below is read directly off that
# paper, not invented.
# ---------------------------------------------------------------------------

PAPER_ELEMENT_REMOVAL_CS = {
    # symbol: (Sigma_R/rho paper value [cm^2/g], tolerance)
    "C": 0.053,
    "H": 0.205,
    "B": 0.058,
    "W": 0.010,
    "O": 0.044,
    "Fe": 0.020,
    "Al": 0.032,
    "Si": 0.031,
}


@pytest.mark.parametrize("symbol,expected", PAPER_ELEMENT_REMOVAL_CS.items())
def test_element_removal_cross_section_matches_paper(symbol, expected):
    got = element_removal_cross_section(symbol)
    assert got == pytest.approx(expected, abs=0.001)


def test_paraffin_wax_total_removal_cross_section_matches_paper():
    # Table 1: pure paraffin wax (CH2)_n -- C fraction 0.852, H fraction
    # 0.147 (the paper's own printed fractions, not renormalized), density
    # such that partial densities are C=0.809, H=0.140 g/cm^3 -- i.e. bulk
    # density 0.809/0.852 = 0.94955 g/cm^3. Paper's own total Sigma_R: 0.071.
    density = 0.809 / 0.852
    composition = {"C": 0.852, "H": 0.147}
    sigma_r = removal_cross_section(composition, density)
    assert sigma_r == pytest.approx(0.071, abs=0.002)


def test_neutron_hvl_and_mfp_formulas():
    sigma_r = 0.1
    assert neutron_hvl(sigma_r) == pytest.approx(6.93)
    assert neutron_mfp(sigma_r) == pytest.approx(10.0)


def test_hydrogen_has_far_higher_mass_removal_cross_section_than_lead():
    # Sanity check independent of exact sourced values: hydrogen has by far
    # the highest MASS (per-gram) removal cross section of any element --
    # the standard, well-known qualitative fact that makes hydrogenous
    # materials good fast-neutron moderators/shields pound-for-pound, unlike
    # dense high-Z metals such as lead, which are excellent gamma shields
    # but poor fast-neutron ones on a mass basis. (The MACROSCOPIC,
    # density-weighted Sigma_R for solid lead metal can still come out
    # numerically larger than water's per unit thickness, since lead is
    # ~11x denser -- that's a real, separate effect, not tested here; this
    # test checks the per-gram comparison, which is the standard one in the
    # removal-cross-section literature.)
    assert element_removal_cross_section("H") > 10 * element_removal_cross_section("Pb")


# ---------------------------------------------------------------------------
# Effective conductivity (Ceff): no local EpiXS ground truth was available
# for this quantity (unlike Zeff/Neff/Zeq), so these are formula/sanity
# checks against the cited literature equations themselves, not an
# independent validation.
# ---------------------------------------------------------------------------


def test_effective_conductivity_is_positive_and_scales_with_density():
    from zeff import effective_conductivity

    water = {"H": 0.111898, "O": 0.888102}
    c1 = effective_conductivity(water, density_g_cm3=1.0, energy_mev=1.0)
    c2 = effective_conductivity(water, density_g_cm3=2.0, energy_mev=1.0)
    assert c1 > 0
    assert c2 == pytest.approx(2 * c1, rel=1e-9)  # Ceff is linear in rho by construction
