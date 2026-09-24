import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "epixs"))

import pytest

from macro import MacroError, run_macro_file


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_macro_mixture_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    macro = _write(tmp_path, "mix.mac", """
        /material/name TESTMIX
        /material/density 1.33
        /material/addElement C 0.5166
        /material/addElement H 0.0745
        /material/addElement O 0.3281
        /material/addElement N 0.0808

        /energies 1.0
        /depths 1,5

        /output/massAttenuation mac_out_mass.csv
        /output/ebf mac_out_ebf.csv
        /output/neutron mac_out_neutron.csv
    """)
    run_macro_file(macro)

    with open(tmp_path / "mac_out_mass.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    # Cross-checked in this project's session history against CHON's own
    # recorded XCOM value for this exact composition: 0.06832 cm^2/g.
    assert float(rows[0]["Total w/ coherent (cm2/g)"]) == pytest.approx(0.06832, abs=1e-4)

    assert (tmp_path / "mac_out_ebf.csv").exists()
    assert (tmp_path / "mac_out_neutron.csv").exists()


def test_macro_single_element_formula(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    macro = _write(tmp_path, "pb.mac", """
        /material/formula Pb
        /energies 1.0
        /output/massAttenuation pb_out.csv
    """)
    run_macro_file(macro)
    with open(tmp_path / "pb_out.csv") as f:
        rows = list(csv.DictReader(f))
    assert float(rows[0]["Total w/ coherent (cm2/g)"]) == pytest.approx(0.07102, abs=1e-4)


def test_macro_rejects_formula_and_addelement_together(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    macro = _write(tmp_path, "bad.mac", """
        /material/formula Pb
        /material/addElement O 0.5
        /energies 1.0
        /output/massAttenuation x.csv
    """)
    with pytest.raises(MacroError):
        run_macro_file(macro)


def test_macro_missing_prerequisite_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    macro = _write(tmp_path, "bad2.mac", """
        /material/formula Pb
        /output/massAttenuation x.csv
    """)
    with pytest.raises(MacroError, match="needs /energies"):
        run_macro_file(macro)


def test_macro_unknown_directive_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    macro = _write(tmp_path, "bad3.mac", "/nonsense/thing 5\n")
    with pytest.raises(MacroError, match="unknown directive"):
        run_macro_file(macro)
