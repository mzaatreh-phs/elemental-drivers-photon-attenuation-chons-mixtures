# XCOM / G-P results for the exact CHONS compositions (epixs_cli), 19 Sep 2026

Script `analysis/xcom_set.py`; data in this folder (`zeff.csv`, `process_fractions.csv`, `crossover.csv`, `exposure_buildup.csv`, `absorption_buildup.csv`). "old" = the original CHON-only Table 2 compositions, "new" = exact compositions including sulfur. Validation: epixs_cli reproduces the EpiXS exports (Zeff median 0.04% / max 0.08%; EBF/EABF median 0.013% / max 0.78%). With the old compositions epixs_cli reproduces the manuscript's EBF numbers (e.g. 0.1 MeV/40 mfp Mix 5 94,131 vs 93,955; 15 MeV Mix 10 7.9416 vs 7.941).

## Effective atomic number (direct method)
| E (MeV) | new range | old range | water | PMMA |
|---|---|---|---|---|
| 0.01 | 6.987 (Mix 10) - 7.394 (Mix 4) | 6.587 - 6.986 | 7.571 | 6.633 |
| 0.03 | 4.586 (Mix 2) - 4.961 (Mix 4) | 4.389 - 4.772 | 4.586 | 4.365 |
| 0.1 | 3.696 (Mix 2) - 3.988 (Mix 4) | 3.669 - 3.949 | 3.403 | 3.643 |
| 1 | 3.646 (Mix 2) - 3.932 (Mix 4) | 3.627 - 3.904 | 3.335 | 3.601 |
| 15 | 4.115 (Mix 2) - 4.418 (Mix 4) | 4.082 - 4.380 | 3.951 | 4.054 |

- Mixture 10 (no sulfur) now has the LOWEST Zeff at 0.01 MeV; the BSA/sulfur-rich mixtures are highest.
- Water is above the CHONS band at 0.01 MeV, within at 0.02, below from 0.03 MeV. PMMA is BELOW the CHONS band at every energy checked (it was inside the old band at 0.01 MeV).

## Photoelectric-incoherent crossover
new 25.34 - 26.18 keV (Mix 1: 25.73); old 24.01 - 25.43 keV (Mix 1: 24.56).

## Process fractions, Mix 1 (new | old)
0.01 MeV: PE 91.4% / INC 3.8% / COH 4.8%  |  90.5 / 4.3 / 5.3.  0.1 MeV: PE 1.3%, INC 96.0%  |  1.0, 96.3.  15 MeV: pair 33.1%  |  32.9.  Coherent peak 12.1% at 0.03 MeV.

## XCOM MAC inter-mixture spread (theoretical)
| E (MeV) | new | old |
|---|---|---|
| 0.01 | 8.36% (Mix 3 max, Mix 2 min) | 17.12% (Mix 10 / Mix 5) |
| 0.03 | 3.08% (Mix 4 max, Mix 10 min) | 6.05% (Mix 10 / Mix 5) |
| 0.1 | 0.88% | 0.93% |
| 1 | 1.12% | 1.10% |
| 10 | 0.39% | 0.65% |
| 15 | 0.48% | 1.01% |

## EBF at 40 mfp (G-P)
| E (MeV) | new | old |
|---|---|---|
| 0.04 | max Mix 10 839.7, min Mix 4 531.9, spread 57.9% | max Mix 5 2194, min Mix 10 840, spread 161% |
| 0.1 | max Mix 10 49,994, min Mix 4 36,164, spread 38.2% | max Mix 5 94,131, min Mix 10 50,014, spread 88.2% |
| 15 | Mix 10 7.942 vs Mix 2 7.871, spread 0.9% | 7.942 vs 7.858, 1.1% |

## What this means for the paper
1. The low-energy story changes: with sulfur, oxygen is no longer the controlling element for MAC at 0.01-0.03 MeV (Mix 10, the most oxygen-rich, is now the LOWEST MAC at 0.03 MeV). Sulfur (Z=16) and BSA content matter. The MAC spread at 0.01 MeV halves (17.1% -> 8.4%).
2. The buildup ranking reverses: Mix 10 now has the HIGHEST EBF at 0.04 and 0.1 MeV (it was the lowest), and the spread at 0.1 MeV/40 mfp drops from 88% to 38%.
3. Zeff conclusions about PMMA and water change (see above); the crossover moves up by about 1 keV.
4. The "exclusively CHON, Z = 1-8, oxygen governs" framing has to be rewritten around CHONS.
These are theoretical XCOM/G-P results; the Geant4 MAC data (running) must confirm the MAC part.
