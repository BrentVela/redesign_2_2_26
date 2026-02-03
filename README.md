# ULTIMATE Alloy Redesign

This repo captures reproducible analysis artifacts for Dr. Karaman's ULTIMATE alloy redesign work.

## Status
- Inputs currently present: `1_find_stoich_props_parallel_for_Class.py`, `TC_Equilibrium_Module.py`, `make_trajectory.py`.
- Missing runtime dependencies: `strength_model` module and `SPACE_v0.par` data file.

## Quick Start
1. Add required inputs:
   - `SPACE_v0.par` (composition space dataset)
   - `strength_model.py` (Curtin model implementation)
2. Run properties extraction:
   - `python 1_find_stoich_props_parallel_for_Class.py`
3. Outputs are written to `CalcFiles/` and CSV files are produced.

## Reproducibility
- Work log: `analysis/LOG.md`
- Property analysis notes: `analysis/properties_analysis.md`

