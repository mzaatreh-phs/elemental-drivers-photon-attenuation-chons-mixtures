# EpiXS-CLI

A command-line (terminal) implementation of the EpiXS algorithm (Hila et al.,
2021) used in this study to compute the XCOM mass attenuation coefficients,
the effective atomic number (Zeff), the equivalent atomic number (Zeq), and
the geometric-progression (G-P) exposure buildup factors (EBF) of the ten
CHONS mixtures.

For any element, compound, or weight-fraction mixture it computes:

- mass attenuation coefficients: coherent, incoherent, photoelectric, pair
  production, and total, from NIST XCOM v3.1;
- linear attenuation coefficient, HVL, TVL, and MFP;
- Zeff, Zeq, and electron density;
- exposure (EBF) and energy-absorption (EABF) buildup factors by the G-P
  method of ANSI/ANS-6.4.3-1991;
- the fast-neutron removal cross section.

---

## 1. Requirements

- Linux or macOS terminal (on Windows, use WSL)
- Python 3.9 or newer (standard library only)
- `gfortran` to compile NIST XCOM the first time

Install `gfortran` if it is missing:

```bash
sudo apt install gfortran        # Debian / Ubuntu
brew install gcc                 # macOS (Homebrew)
```

## 2. Get the program

```bash
git clone https://github.com/mzaatreh-phs/elemental-drivers-photon-attenuation-chons-mixtures.git
cd elemental-drivers-photon-attenuation-chons-mixtures/epixs_cli
chmod +x run
```

The first time you run it, `./run` compiles XCOM automatically. It prints
`XCOM binary not built yet -- compiling...` and takes a few seconds. To
compile by hand instead:

```bash
cd vendor/xcom
gfortran -w -std=legacy -o xcom XCOM.f
cd ../..
```

## 3. Reproduce the values of the article

Each file `examples/mixNN.mac` contains one mixture with the elemental mass
fractions and density of Table 2 of the article, the 37 energies from 0.015
to 15 MeV, and the depths 5-40 mfp used in the article.

Run one mixture:

```bash
./run examples/mix01.mac
```

Run all ten mixtures:

```bash
for f in examples/mix*.mac; do ./run "$f"; done
```

Each run writes four CSV files into the current folder, for example for
Mixture 1:

| File | Content |
|---|---|
| `mix01_mass_attenuation.csv` | MAC by process and total (cm^2/g) |
| `mix01_hvl_tvl_mfp.csv` | LAC (1/cm), HVL, TVL, MFP (cm) |
| `mix01_zeff_zeq.csv` | Zeff, Zeq, electron density |
| `mix01_ebf.csv` | EBF at 5, 10, ..., 40 mfp |

Check: for Mixture 4 at 0.1 MeV and 40 mfp, `mix04_ebf.csv` gives
EBF = 36164, the value reported in the article.

## 4. Your own material with a macro file (batch mode)

Copy `sample.mac`, edit it, and run it:

```bash
cp sample.mac my_material.mac
nano my_material.mac          # or any text editor
./run my_material.mac
```

A macro file has one directive per line; `#` starts a comment:

```
# a mixture from elemental weight fractions (normalized automatically)
/material/name      MyMixture
/material/density   1.33
/material/addElement C 0.51
/material/addElement H 0.07
/material/addElement O 0.33
/material/addElement N 0.08
/material/addElement S 0.01

# or a single element or compound instead of addElement lines:
#/material/formula  H2O

/energies 0.1,0.5,1,5,10         # MeV
/depths   1,5,10,20,40           # mean free paths, for EBF/EABF

/output/massAttenuation my_mac.csv
/output/hvl             my_hvl.csv
/output/zeff            my_zeff.csv
/output/ebf             my_ebf.csv
/output/eabf            my_eabf.csv
/output/neutron         my_neutron.csv
```

| Directive | Needs |
|---|---|
| `/output/massAttenuation` | `/energies` |
| `/output/hvl` | `/energies`, `/material/density` |
| `/output/zeff` | `/energies`, `/material/density` |
| `/output/ebf`, `/output/eabf` | `/energies`, `/depths` |
| `/output/neutron` | `/material/density` |

## 5. Interactive menu

Run without arguments:

```bash
./run
```

The menu asks what to compute:

```
Main menu -- what would you like to compute?
    1. Mass attenuation coefficient
    2. Effective atomic number (Zeff), electron density & conductivity
    3. Exposure buildup factor (EBF)
    4. Energy absorption buildup factor (EABF)
    5. Fast neutron removal cross section
    6. Exit
```

It then asks for the material, either a symbol or formula such as `Pb` or
`H2O`, or a mixture entered component by component with weight fractions,
followed by the energies in MeV and, for buildup factors, the depths in mfp.
Results are printed in the terminal and can be saved to a CSV file.

## 6. Tests

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install pytest
python -m pytest tests/ -q
```

On Debian or Ubuntu, `python3 -m venv` needs the `python3-venv` package
(`sudo apt install python3-venv`). The current version passes all 859 tests.

The tests compare the buildup factors with reference values exported from
the original EpiXS program, which are embedded in the test files.

## 7. Layout

| Path | Content |
|---|---|
| `run` | launcher: compiles XCOM if needed, then starts the menu or a macro |
| `epixs/xcom_engine.py` | Python driver for the compiled XCOM program |
| `epixs/gp_data.py`, `epixs/gp_calculator.py` | G-P coefficients and buildup formula |
| `epixs/zeff.py` | Zeff, Zeq, electron density |
| `epixs/materials.py` | formula parsing and mixture composition |
| `epixs/macro.py` | macro-file batch mode |
| `epixs/cli.py` | interactive menu |
| `examples/` | the ten CHONS mixtures of the article, and lead examples |
| `vendor/xcom/` | NIST XCOM v3.1 Fortran source and data |

## 8. Sources and credits

- Algorithm: F. C. Hila et al., EpiXS: A Windows-based program for photon
  attenuation, dosimetry and shielding based on EPICS2017 (ENDF/B-VIII) and
  EPDL97 (ENDF/B-VI.8), Radiation Physics and Chemistry 182 (2021) 109331.
- Photon cross sections: M. J. Berger et al., XCOM: Photon Cross Sections
  Database, NIST Standard Reference Database 8 (XGAM); XCOM v3.1 Fortran
  program, distributed by NIST.
- Buildup coefficients: ANSI/ANS-6.4.3-1991, Gamma-Ray Attenuation
  Coefficients and Buildup Factors for Engineering Materials, American
  Nuclear Society.
