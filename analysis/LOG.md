# Work Log

All timestamps in UTC.

## 2026-02-03T16:49:11Z
- Sanitized filenames and removed Zone.Identifier artifacts.
- Initialized reproducibility scaffolding: `.gitignore`, `README.md`, `analysis/LOG.md`, `analysis/properties_analysis.md`.

- Initialized git repository and committed baseline files on branch `main` (commit aa99560).
## 2026-02-03T16:59:20Z
- Ran `make_trajectory.py` to generate `trajectory.csv` and `SPACE_v0.par`.
- Ran `1_find_stoich_props_parallel_for_Class.py` with `strength_model.py` to compute properties.
- Output written to `CalcFiles/STOIC_OUT_0.csv`.
## 2026-02-03T17:06:58Z
- Added `plot_properties_vs_trajectory.py` (style matched to `plot_ys_vs_temperature.py`).
- Generated `properties_vs_trajectory.pdf` from `CalcFiles/STOIC_OUT_0.csv` (with `--no-tex`).
- Removed Zone.Identifier artifacts for plotting and strength model files.
## 2026-02-03T17:11:20Z
- Updated default plotted properties in `plot_properties_vs_trajectory.py` (removed entropy/hardness/Zener/bulk/shear; added Pugh ratio).
## 2026-02-03T17:15:07Z
- Added output format and layout flags to `plot_properties_vs_trajectory.py` to speed LaTeX rendering and allow SVG/PNG.
## 2026-02-03T17:20:23Z
- Added one-per-figure mode and TeX prewarm option to `plot_properties_vs_trajectory.py`.
## 2026-02-03T17:24:59Z
- Disabled TeX by default for trajectory plotting and switched labels to Unicode subscripts.
## 2026-02-03T17:26:08Z
- Generated one-per-figure trajectory property plots (no TeX) in `plots/`.
## 2026-02-03T17:29:18Z
- Updated YS labels to use Unicode subscripts with degree symbol.
- Generated 1x6 stacked plot: `properties_vs_trajectory_1x6.pdf`.
## 2026-02-03T17:36:04Z
- Updated default plots to use YS at test temperature (1300C) and added reference lines for density/YS.
- Regenerated `properties_vs_trajectory_1x6.pdf`.
## 2026-02-03T17:47:45Z
- Merged EQUIL outputs with stoichiometric properties into `CalcFiles/STOIC_EQUIL_OUT.csv` (11 rows).
## 2026-02-03T18:41:43Z
- Added in-repo rerun script for missing Thermo-Calc equilibria: `rerun_equil_missing.py`.
## 2026-02-03T18:48:18Z
- Added TC25B_HOME auto-detection to TC scripts (defaults to /opt/Thermo-Calc/2025b).
## 2026-02-03T19:07:30Z
- Stitched missing equilibrium rows and regenerated `CalcFiles/STOIC_EQUIL_OUT.csv`.
- Plotted 500C BCC fraction sum vs V: `bcc_fraction_vs_v.pdf`.
## 2026-02-04
- Set 3_creep_resistance.py grain size to paper value (1.2 mm); docstring retained.
- Regenerated parity data/plot for TCHEA7+MOBHEA3, rate-limiting DI: fig8_parity_rate_limit_TCHEA7_MOBHEA3.pdf.
- Created shifted parity plots using systematic offset (+0.386 log10):
  - fig8_parity_rate_limit_TCHEA7_MOBHEA3_shifted.pdf
  - fig8_parity_rate_limit_TCHEA7_MOBHEA3_100h_pct_shifted.pdf
- Predicted shifted 1200C, 150 MPa, 100 h creep strain along trajectory:
  - CalcFiles/creep_strain_1200C_150MPa_100h_VG_shifted.csv
  - creep_strain_1200C_150MPa_100h_vs_v_pct_shifted.pdf
  - creep_strain_1200C_150MPa_100h_vs_v_pct_shifted_no_max.pdf (max removed)
- Added composite a-priori creep plot script (min of NH, Coble, PL5): plot_creep_composite_vs_v.py
  - CalcFiles/creep_composite_1200C_vs_v.csv
  - creep_composite_1200C_vs_v.pdf
  - creep_composite_1200C_vs_v_100h_pct.pdf
- Generated NH-only 100 h % plot: creep_NH_1200C_vs_v_100h_pct.pdf
