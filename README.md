# Data and analysis for: Lead-Free CHONS-Based Biomolecular Mixtures for Photon Attenuation

Supporting data and code for the article *Lead-Free CHONS-Based Biomolecular Mixtures for Photon Attenuation: A Geant4 Monte Carlo and Pearson Correlation Investigation* (M. Y. AlZaatreh and N. Z. Noor Azman). Ten idealized mixtures of alpha-cellulose, bovine serum albumin (BSA), sorbitol and stearic acid (C, H, O, N and 0-1.5 wt% S) were simulated with Geant4 over 0.01-15 MeV.

**License:** code (analysis scripts and `epixs_cli/`) under the MIT License (`LICENSE`); data, macros and figure data under CC BY 4.0 (`LICENSE-DATA`). Please cite the article when using this material.

## What is here
| Path | Content | Raw or derived |
|---|---|---|
| `raw/<job>/`, `macros/<job>.mac` | Main campaign: 760 Geant4 jobs (10 mixtures x 38 energies x 2 seeds x 10^6 photons): result CSV, console log, exact macro | raw |
| `raw_exchange/<job>/`, `macros_exchange/` | Elemental-exchange campaign: 96 Geant4 jobs (base composition and five variants x 8 energies x 2 seeds x 10^6 photons) | raw |
| `inputs/` | Exact elemental compositions, densities and energy grids used to build the macros (`inputs_exact.json`, `densities_fixed.json`) and the exchange design (`exchange_design.json`) | input |
| `analysis/` | Python/shell scripts that pool the seeds, compute uncertainties, XCOM/G-P quantities, tables, sensitivity analysis and figures, together with their derived result tables (`*/results`, `xcom/`, `sensitivity/`, `broadbeam/`) | scripts + derived |
| `figure_data/` | The numbers behind every figure of the article, one folder per figure (see its README), and the Excel workbook of all sheets | derived |
| `epixs_cli/` | EpiXS-CLI, the authors' terminal implementation of the EpiXS algorithm used for the XCOM coefficients, Z_eff, Z_eq and G-P buildup factors, with the ten mixtures of the article as ready-to-run examples (see `epixs_cli/README.md`) | software |
| `SHA256SUMS.txt` | SHA-256 checksum of every file in this repository | integrity |

Everything under `raw/` and `raw_exchange/` is unmodified simulation output. Pooled values, tables and figure data are derived and can be regenerated from the raw counts with the scripts in `analysis/`.

## Which script produces what
| Article item | Script (in `analysis/`) | Output |
|---|---|---|
| Pooled MAC, uncertainties, seed check, and validation against XCOM | `g4/analyze.py`, `g4/g4analysis.py`, `g4/verify_independent.py` (pure-Python re-derivation) | `g4/results/` |
| Inter-mixture spread and Pearson correlation tables | `g4/make_tables.py` | `g4/results/table4_new.tex`, `table5_new.tex` |
| Zeff, process fractions, photoelectric-incoherent crossover, and G-P exposure buildup | `xcom_set.py` | `xcom/` |
| Carbon-referenced elemental contrasts, reduced-collinearity ensemble, and range test | `sensitivity/elemental_sensitivity.py`, `sensitivity/range_test.py` | `sensitivity/` |
| Geant4 elemental-exchange confirmation | `sensitivity/exchange_confirmation.py`, `sensitivity/gen_exchange_jobs.py` | `sensitivity/exchange_confirmation.csv` |
| Equal-thickness and equal-areal-mass broad-beam response | `broadbeam/equal_thickness.py` | `broadbeam/` |
| Moisture sensitivity with 10 wt% water added to every mixture (supplementary analysis, not reported in the article) | `moisture/moisture_sensitivity.py`, `moisture/moisture_patch.py` | `moisture/` |
| Elemental contrast and validation-residual figures | `origin_style/make_new_figures.py` | `figure_data/fig_elemental_contrast/`, `figure_data/fig_validation_residuals/` |
| Composition scenarios, residual correlations, and dominant-element tables | `sensitivity/scenario_correlations.py`, `sensitivity/scenario_tables_patch.py` | `sensitivity/scenario_correlations.csv` |
| Scoping survey of low-Z literature used to position the article | `literature_survey/survey.py` | `literature_survey/LITERATURE_SURVEY_2026-09-21.md` (approximate keyword statistics, not a systematic review) |
| Plotted data of every figure | `g4/export_figure_data.py` | `figure_data/` |
| Figures in the style of the article | `origin_style/make_core_figures.py`, `origin_style/make_other_figures.py` | (images not included) |
| Numbers quoted in the Results text | `results_text/numbers_*.py` | printed values |

## Reproducing the analysis
- The scripts contain absolute paths of the authors' workstation (`~/CHON/recalc_2026-09-19/...`); adjust the path constants at the top of each script (`R`, `A`, `DEST`, `FD`) before running.
- Python 3 with numpy, scipy, matplotlib and openpyxl.
- `xcom_set.py` and `broadbeam/equal_thickness.py` need **EpiXS-CLI**, the authors' command-line implementation of the EpiXS algorithm (compiled NIST XCOM v3.1 Fortran cross sections plus ANSI/ANS-6.4.3-1991 G-P coefficients). It is included in `epixs_cli/`; see `epixs_cli/README.md` for terminal usage. Point the `sys.path` lines at the top of these two scripts to that folder. Its XCOM- and G-P-derived outputs are also included in `analysis/xcom/` and `analysis/broadbeam/`.
- The Geant4 executable was built from the authors' local `brks` project (source snapshot identified below); the macros in this repository fully define every run.

## Elemental-exchange campaign (`raw_exchange/`)
Job naming `<seedset>_<variant>_E<idx>`: seed set `s1`/`s2`; variant `base` (mean composition of the ten mixtures, density 1.352033 g/cm3), `dO` (+5 wt% O replacing C), `dN` (+5 wt% N), `dS` (+1.5 wt% S), `dH` (+1.5 wt% H) or `dALL` (all four together); energy index `00`-`07` = 0.01, 0.02, 0.03, 0.05, 0.1, 0.3, 1, 10 MeV. Setup and counting are identical to the main campaign. The exact compositions are in `inputs/exchange_design.json` and in each macro.

## Main campaign: job naming
`<seedset>_Mix_<m>_E<idx>`, for example `s1_Mix_5_E08`: seed set `s1` or `s2` (two independent random-number streams), mixture number `m` (1-10), energy index `idx` (table below). 10 mixtures x 38 energies x 2 seed sets = **760 jobs**, all completed (0 failed).

## Simulation setup (identical for all jobs)
- Geant4 11.4.2 (`geant4-11-04-patch-02`), physics list `G4EmStandardPhysics_option4+G4DecayPhysics+G4RadioactiveDecayPhysics`, production cut 0.001 mm.
- Primary photons per job: **1,000,000**, monoenergetic, one energy per job, one thread per job.
- Narrow-beam geometry: point source, lead collimator with a 0.3 cm hole, sample slab **1.0 cm thick (fixed for every energy and mixture)**, 10 cm source-to-sample and sample-to-detector gaps, collimator blocks 22.7 cm in diameter.
- `n_incident` counts primary photons entering the slab in the forward direction (once per event). `n_transmitted_strict` counts those primary photons leaving the far face with direction cosine > 0.999 (about 2.6 degrees). `n_transmitted_forward` counts primary photons leaving the far face in the forward hemisphere. Counts refer to the primary photon track only.
- Random seeds: `/random/setSeeds` in each macro (values differ for every job).
- Executable SHA-256: `a1d52dc417f3be5fab09028aa8e8a5610118c100b80b5b26c0c9878cdfe8839a`. Source snapshot: base commit `5feccf955c97488f7a6a7e62054ecdd371330ef4` of the authors' local `brks` project plus uncommitted local modifications; the SHA-256 of the list of source-file hashes is `e0236c4273a363ef1dddd2228499e140dff9a437a40c092093dc21d948c61249`. (The code itself is not part of this repository.)

## Derived quantity used by the authors
`T = n_transmitted_strict / n_incident`, `MAC = -ln(T) / (rho * x)` with `x = 1.0 cm` and `rho` the density below; binomial standard uncertainty `u(MAC) = sqrt((1-T)/(n_incident*T)) / (rho*x)`. The other columns of the CSV are computed inside the simulation code from the same counts; the `Zeff` and `IdealMAC` columns there are the code's own values and are not used by the authors' analysis.

## Mixtures (material definitions exactly as in the macros)
Constituent proportions: alpha-cellulose / BSA / sorbitol / stearic acid (wt%). Elemental mass fractions were calculated from the chemical formulas of the anhydrous constituents; BSA = C2934 H4581 N781 O897 S39 from the UniProt P02769 mature chain (17 disulfide bonds). Density by the additive-volume rule from 1.50 (cellulose), 1.33 (BSA), 1.49 (sorbitol), 0.94 (stearic acid) g/cm3.

| Mix | cellulose / BSA / sorbitol / stearic (wt%) | C | H | O | N | S | Density (g/cm3) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 30 / 50 / 10 / 10 | 0.514264 | 0.073925 | 0.320019 | 0.082377 | 0.009415 | 1.334335 |
| 2 | 30 / 50 / 0 / 20 | 0.550701 | 0.078934 | 0.278573 | 0.082377 | 0.009415 | 1.267900 |
| 3 | 30 / 50 / 20 / 0 | 0.477826 | 0.068916 | 0.361466 | 0.082377 | 0.009415 | 1.408117 |
| 4 | 37.5 / 62.5 / 0 / 0 | 0.498385 | 0.066778 | 0.320097 | 0.102971 | 0.011769 | 1.389034 |
| 5 | 0 / 80 / 10 / 10 | 0.540145 | 0.076138 | 0.236850 | 0.131803 | 0.015065 | 1.290321 |
| 6 | 20 / 60 / 10 / 10 | 0.522891 | 0.074662 | 0.292296 | 0.098852 | 0.011298 | 1.319334 |
| 7 | 40 / 40 / 10 / 10 | 0.505636 | 0.073187 | 0.347743 | 0.065902 | 0.007532 | 1.349681 |
| 8 | 50 / 30 / 10 / 10 | 0.497009 | 0.072450 | 0.375466 | 0.049426 | 0.005649 | 1.365385 |
| 9 | 60 / 20 / 10 / 10 | 0.488382 | 0.071712 | 0.403189 | 0.032951 | 0.003766 | 1.381458 |
| 10 | 80 / 0 / 10 / 10 | 0.471128 | 0.070237 | 0.458635 | 0.000000 | 0.000000 | 1.414767 |

## Energy grid (MeV)
| idx | E (MeV) | role |
|---:|---:|---|
| 00 | 0.01 | production |
| 01 | 0.015 | production |
| 02 | 0.02 | production |
| 03 | 0.03 | production |
| 04 | 0.04 | production |
| 05 | 0.05 | production |
| 06 | 0.06 | production |
| 07 | 0.08 | production |
| 08 | 0.1 | production |
| 09 | 0.15 | production |
| 10 | 0.2 | production |
| 11 | 0.3 | production |
| 12 | 0.356 | validation |
| 13 | 0.4 | production |
| 14 | 0.5 | production |
| 15 | 0.511 | validation |
| 16 | 0.6 | production |
| 17 | 0.662 | validation |
| 18 | 0.8 | production |
| 19 | 1 | production |
| 20 | 1.173 | validation |
| 21 | 1.33 | validation |
| 22 | 1.5 | production |
| 23 | 2 | production |
| 24 | 2.51 | validation |
| 25 | 3 | production |
| 26 | 4 | production |
| 27 | 5 | production |
| 28 | 6 | production |
| 29 | 7 | production |
| 30 | 8 | production |
| 31 | 9 | production |
| 32 | 10 | production |
| 33 | 11 | production |
| 34 | 12 | production |
| 35 | 13 | production |
| 36 | 14 | production |
| 37 | 15 | production |

## CSV columns
`material_id, density_g_cm3, Energy(MeV), n_incident, n_transmitted_strict, n_transmitted_forward, T_strict, T_strict_uncertainty, MAC_strict(cm2/g), T_forward, MAC_forward(cm2/g), IdealMAC(cm2/g), LAC(1/cm), HVL(cm), TVL(cm), MFP(cm), Zeff, RPE(%), OpticalDensity, CrossSection(barn), Thikness(cm), geant4_version, physics_list, seed, run_id, timestamp`
(`Thikness(cm)` is spelled as written by the code.) `density_g_cm3` is printed by the code with reduced precision; the density actually used is the value in the job's macro.
