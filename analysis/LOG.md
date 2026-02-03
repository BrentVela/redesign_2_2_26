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
