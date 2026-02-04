# ULTIMATE Alloy Redesign

This repo captures reproducible analysis artifacts for Dr. Karaman's ULTIMATE alloy redesign work.

## Status
- Inputs currently present: `1_find_stoich_props_parallel_for_Class.py`, `TC_Equilibrium_Module.py`, `make_trajectory.py`.
- Missing runtime dependencies: `strength_model` module and `SPACE_v0.par` data file.

## Quick Start
1. Add required inputs: `SPACE_v0.par` (composition space dataset), `strength_model.py` (Curtin model implementation).
2. Run properties extraction: `python 1_find_stoich_props_parallel_for_Class.py`.
3. Outputs are written to `CalcFiles/` and CSV files are produced.

## Thermo-Calc / tcpython Notes
- A working Thermo-Calc install is required for `tc_python`.
- The conda env used here is `tcpython` and needs `PYTHONPATH` set to the repo root.
- `TC25B_HOME` should point to the local Thermo-Calc install (e.g., `/home/vela/thermocalc/2025b`).
- For the TC property run used here, we invoked `analysis/run_tc_properties.py` via `conda run -n tcpython --cwd /home/vela/projects/ult/redesign_2_2_26 python analysis/run_tc_properties.py`.
- In the Codex sandbox, Thermo-Calc startup can hang; running outside the sandbox was required.

## Reproducibility
- Work log: `analysis/LOG.md`
- Property analysis notes: `analysis/properties_analysis.md`
