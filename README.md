# CHON shielding: raw Geant4 narrow-beam attenuation data

Raw output of the Geant4 campaign for ten CHONS biomolecular mixtures (alpha-cellulose, bovine serum albumin, sorbitol, stearic acid).
Produced on 2026-09-19 (16:34-21:47 UTC+03:00). **Private repository, pre-publication.** License and citation to be added by the authors before this repository is made public.

## What is here
| Path | Content |
|---|---|
| `raw/<job>/MAC_results_Mix_<m>.csv` | one-row result file written by the simulation for that job |
| `raw/<job>/run.log` | full Geant4 console log of that job |
| `macros/<job>.mac` | the exact macro that was run (material, density, thickness, energy, seed) |
| `SHA256SUMS.txt` | SHA-256 checksums of every file above |

Nothing in this repository is derived: no pooled values, tables or figures. Derived quantities are computed from the counts as described below.

## Job naming
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
