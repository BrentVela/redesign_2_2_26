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
