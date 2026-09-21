figs_new/fig_composition_sensitivity.png - two panels, x log axis (Energy_MeV).
PANEL A (left, log y): MAC inter-mixture spread (%) -> panelA_MAC_spread_percent.csv.
   Solid line with markers = Geant4_spread; shaded band between Geant4_noise_band_16pct and _84pct (68% Monte Carlo noise band); dashed black = XCOM_spread.
PANEL B (right, y from -1.05 to 1.05): Pearson r between MAC and the mass fraction of the element -> panelB_Pearson_r_with_MAC.csv.
   Solid + markers = *_Geant4 columns, dashed same colour = *_XCOM columns (Oxygen, Sulfur, Hydrogen).
Note: sulfur, nitrogen and BSA fraction are perfectly collinear in this ten-mixture library, so the sulfur curve is also the nitrogen and BSA curve.
Above 1 MeV the Geant4 spread is noise-limited; use the XCOM columns there.
