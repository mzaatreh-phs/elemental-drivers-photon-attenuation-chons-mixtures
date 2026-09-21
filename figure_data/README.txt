FIGURE DATA - recalculated CHON study (Geant4 + XCOM/G-P), regenerate with ~/CHON/recalc_2026-09-19/analysis/g4/export_figure_data.py
LAYOUT: one folder per figure. Folder name = PNG name in ../figs_new/ (e.g. fig_mac_lac/ <-> figs_new/fig_mac_lac.png).
Inside: panelA_..., panelB_... = panels in reading order (left->right, top->bottom); each folder has its own README.txt with the x column, y columns,
axis scales and what is a curve / error bar / shaded band. All files are CSV with header row; "wide" = Energy_MeV then one column per curve.
figure_data_all_sheets.xlsx = same tables, sheet 'INDEX' lists which sheet is which file. Mixture composition: Table2_mixture_composition_and_density.csv.
Not in any figure: fig1_geometry (schematic, no data).
Uncertainties are binomial from the pooled counts (2 x 10^6 photons per point); above 1 MeV the Geant4 composition spread is noise-limited (use XCOM there).

FIGURE FOLDER                              SOURCE                        CONTENT

fig_mac_lac/                               Geant4                        (see fig_mac_lac/README.txt)
   panelA_MAC_curves.csv                                             MAC (Geant4) - the curves of the left panel
   panelA_MAC_errorbars_1sigma.csv                                   MAC 1-sigma binomial uncertainty (error bars of left panel; smaller than the markers in the PNG)
   panelA_MAC_XCOM_reference.csv                                     XCOM total MAC with coherent scattering, same energies - NOT plotted in the PNG, optional overlay
   panelB_LAC_curves.csv                                             LAC = MAC x density (Geant4) - the curves of the right panel
   panelB_LAC_errorbars_1sigma.csv                                   LAC 1-sigma uncertainty (error bars of right panel)

fig_mfp_hvl_tvl/                           Geant4                        (see fig_mfp_hvl_tvl/README.txt)
   panelA_MFP_cm.csv                                                 mean free path = 1/LAC (Geant4)
   panelB_HVL_cm.csv                                                 half-value layer = ln2/LAC (Geant4)
   panelC_TVL_cm.csv                                                 tenth-value layer = ln10/LAC (Geant4)

fig_rpe/                                   Geant4 (1 cm)                 (see fig_rpe/README.txt)
   curves_attenuation_1cm_percent.csv                                primary-beam attenuation (1 - T) x 100 for the 1.0 cm slab (Geant4) - the curves
   errorbars_1sigma_percent.csv                                      1-sigma of the attenuation (error bars)

fig_validation_deviation/                  Geant4 vs XCOM                (see fig_validation_deviation/README.txt)
   points_deviation_percent.csv                                      signed deviation 100 x (MAC_Geant4/MAC_XCOM - 1) - the plotted points
   errorbars_1sigma_percent.csv                                      error bar length (1 sigma) of each point
   validation_energies_MeV.csv                                       the six validation energies (dotted vertical lines in the PNG)

fig_seed_check/                            Geant4 (2 seeds)              (see fig_seed_check/README.txt)
   z_values_all_points.csv                                           z = (MAC_seed1 - MAC_seed2)/sqrt(u1^2+u2^2) for every mixture and energy - raw data of the histogram
   histogram_and_N01_curve.csv                                       histogram bins (bar height = density) and the standard-normal curve (black line in the PNG)

fig_zeff/                                  XCOM                          (see fig_zeff/README.txt)
   curves_Zeff.csv                                                   Zeff of each mixture (XCOM direct method)

fig_zeff_baseline/                         XCOM                          (see fig_zeff_baseline/README.txt)
   band_and_reference_Zeff.csv                                       shaded band (min..max Zeff over the ten mixtures) + water and PMMA lines

fig_process_fractions/                     XCOM (Mix1)                   (see fig_process_fractions/README.txt)
   Mix1_process_fractions_percent.csv                                Mixture 1: percentage of the total MAC from each process (4 curves)
   crossover_keV_all_mixtures.csv                                    photoelectric = incoherent energy per mixture; the dashed vertical line in the PNG is the Mix1 value

fig_ebf/                                   G-P                           (see fig_ebf/README.txt)
   panel1_5mfp_EBF.csv                                               exposure buildup factor at 5 mfp (G-P), panel 1 of 8
   panel2_10mfp_EBF.csv                                              exposure buildup factor at 10 mfp (G-P), panel 2 of 8
   panel3_15mfp_EBF.csv                                              exposure buildup factor at 15 mfp (G-P), panel 3 of 8
   panel4_20mfp_EBF.csv                                              exposure buildup factor at 20 mfp (G-P), panel 4 of 8
   panel5_25mfp_EBF.csv                                              exposure buildup factor at 25 mfp (G-P), panel 5 of 8
   panel6_30mfp_EBF.csv                                              exposure buildup factor at 30 mfp (G-P), panel 6 of 8
   panel7_35mfp_EBF.csv                                              exposure buildup factor at 35 mfp (G-P), panel 7 of 8
   panel8_40mfp_EBF.csv                                              exposure buildup factor at 40 mfp (G-P), panel 8 of 8

fig_composition_sensitivity/               Geant4 + XCOM                 (see fig_composition_sensitivity/README.txt)
   panelA_MAC_spread_percent.csv                                     inter-mixture MAC spread; line+markers = Geant4_spread, shaded band = 16..84 pct columns, dashed black = XCOM_spread
   panelB_Pearson_r_with_MAC.csv                                     Pearson r of MAC with O, S, H mass fraction; solid = Geant4 (markers), dashed = XCOM

fig_elemental_sensitivity/                 XCOM mixture rule             (see fig_elemental_sensitivity/README.txt)
   panelA_MAC_change_per_1wtpct_replacing_C_percent.csv              percent change of MAC per +1 wt% of the element replacing carbon (XCOM mixture rule, library-mean composition)
   panelB_share_of_MAC_percent_library_mean.csv                      share of the total MAC contributed by each element (%), library-mean composition

fig_library_vs_ensemble_correlation/       XCOM ensemble                 (see fig_library_vs_ensemble_correlation/README.txt)
   panelA_oxygen_r_with_MAC.csv                                      Pearson r of MAC with oxygen mass fraction: ten-mixture library (red circles) vs decorrelated ensemble (blue squares)
   panelB_sulfur_r_with_MAC.csv                                      Pearson r of MAC with sulfur mass fraction: ten-mixture library (red circles) vs decorrelated ensemble (blue squares)
   panelC_hydrogen_r_with_MAC.csv                                    Pearson r of MAC with hydrogen mass fraction: ten-mixture library (red circles) vs decorrelated ensemble (blue squares)

fig_element_importance/                    XCOM ensemble                 (see fig_element_importance/README.txt)
   importance_percent.csv                                            importance of each element for MAC variability (%) in the decorrelated ensemble (standardized-beta based)

fig_exchange_confirmation/                 Geant4 vs XCOM                (see fig_exchange_confirmation/README.txt)
   panelA_G4_vs_XCOM_MAC_change.csv                                  scatter of Geant4 vs XCOM change in MAC caused by each composition exchange (long format: 5 variants x 8 energies)
   panelB_z_vs_energy.csv                                            z = (Geant4 - XCOM)/uncertainty at each energy, one column per variant (wide)

fig_broadbeam_primary_vs_total/            G-P                           (see fig_broadbeam_primary_vs_total/README.txt)
   panelA_equal10cm_primary_only_relative_to_mean_percent.csv        panel A: equal thickness 10 cm, primary only: transmission of each mixture relative to the mean over the ten mixtures (%)
   panelB_equal20gcm2_primary_only_relative_to_mean_percent.csv      panel B: equal areal mass 20 g/cm2, primary only: transmission of each mixture relative to the mean over the ten mixtures (%)
   panelC_equal10cm_with_buildup_relative_to_mean_percent.csv        panel C: equal thickness 10 cm, with buildup: transmission of each mixture relative to the mean over the ten mixtures (%)
   panelD_equal20gcm2_with_buildup_relative_to_mean_percent.csv      panel D: equal areal mass 20 g/cm2, with buildup: transmission of each mixture relative to the mean over the ten mixtures (%)

fig_broadbeam_spread/                      G-P                           (see fig_broadbeam_spread/README.txt)
   panelA_equal10cm_spread_percent.csv                               panel A: equal thickness 10 cm: best-to-worst mixture spread in transmission (%), primary only vs with buildup
   panelB_equal20gcm2_spread_percent.csv                             panel B: equal areal mass 20 g/cm2: best-to-worst mixture spread in transmission (%), primary only vs with buildup

extra_no_figure/                                                         (see extra_no_figure/README.txt)
   element_MAC_pure_elements.csv                                     MAC of the pure elements C,H,N,O,S (XCOM) - basis of the sensitivity map, not plotted
   range_test_importance.csv                                         Pratt variance shares for the 7 range scenarios (robustness test of element importance), not plotted
   exchange_confirmation_full.csv                                    full exchange-confirmation table incl. additivity rows and absolute MACs
   sensitivity_summary.json                                          summary numbers of the sensitivity analysis
   exchange_design.json                                              compositions used for the exchange runs
   broadbeam_summary_all_conditions.csv                              broad-beam spread/best/worst mixture for every thickness and areal mass
   broadbeam_long.csv                                                broad-beam per mixture, energy and condition (long format)
   broadbeam_wide_B_m20gcm2.csv                                      absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_B_m50gcm2.csv                                      absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_B_t10cm.csv                                        absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_B_t20cm.csv                                        absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_broad_m20gcm2.csv                                absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_broad_m50gcm2.csv                                absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_broad_t10cm.csv                                  absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_broad_t20cm.csv                                  absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_primary_m20gcm2.csv                              absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_primary_m50gcm2.csv                              absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_primary_t10cm.csv                                absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass
   broadbeam_wide_T_primary_t20cm.csv                                absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass

(top level, not a figure)
   Table2_mixture_composition_and_density.csv                        exact elemental mass fractions, density and BSA weight fraction of the ten mixtures (manuscript Table 2); not a figure
