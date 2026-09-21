figs_new/fig_broadbeam_primary_vs_total.png - G-P (epixs_cli) exposure buildup, 2 x 2 panels, x log axis (Energy_MeV), y range -30..30, horizontal line at 0.
y = 100 x (T_mixture / mean over the ten mixtures - 1), % (already computed here, ready to plot). Empty cells = slab thinner than 0.5 mfp (equal-thickness 10 cm: E >= 3 MeV; 20 g/cm2: E >= 7 MeV), outside the G-P validity range 0.5-40 mfp, so no value exists (the PNG leaves them blank).
PANEL A (top-left):     equal thickness 10 cm, primary only          PANEL B (top-right):    equal areal mass 20 g/cm2, primary only
PANEL C (bottom-left):  equal thickness 10 cm, with buildup          PANEL D (bottom-right): equal areal mass 20 g/cm2, with buildup
The absolute transmissions these are computed from are in extra_no_figure/ (broadbeam_wide_T_*, broadbeam_wide_B_*).
CAVEAT: G-P accuracy for these low-Z mixtures is not independently verified (buildup cross-check was dropped); state the ranking claims as G-P-based.
Columns Mix1..Mix10 = the ten mixtures (composition/density in Table2_mixture_composition_and_density.csv).
