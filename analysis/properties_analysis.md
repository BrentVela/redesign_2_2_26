# Properties Analysis (Static)

This analysis is based on the code in `1_find_stoich_props_parallel_for_Class.py` without executing it. Actual results will require `SPACE_v0.par` and `strength_model`.

## Inputs and Assumptions
- Composition space file: `SPACE_v0.par` loaded into a pandas DataFrame.
- Elements list (hard-coded): `['Cr','Mo','Nb','Ta','Ti','V','W']` (sorted).
- Test temperature: hard-coded to `1300` (C) for all rows.
- Strength model dependency: `strength_model.model_Control(elements, comp_i, T=...)`.

## Derived Properties (Per Composition)
The script computes weighted-average (rule of mixtures) properties using hard-coded elemental constants.

- Density Avg [g/cm^3]
- MW Avg [g/mol]
- Tm Avg [K]
- VEC Avg (valence electron concentration)
- Elastic constants (C11, C12, C44) and derived:
  - Cauchy pressure: `C12 - C44`
  - Zener ratio: `2*C44 / (C11 - C12)`
  - Universal anisotropy: `(6/5) * (sqrt(Z) - 1/sqrt(Z))^2`
- Bulk modulus Avg, Shear modulus Avg (from elemental data)
- Electronegativity averages (Pauling, Allen)
- Metallic radius (from `Rm` field)
- Configurational entropy (per composition): `Sconf = -sum(x_i * ln(x_i))`

## Strength/Hardness Estimates
Using the Curtin model (`strength_model.model_Control`), the script calculates:
- `Tau_y 0`, `delta Eb`, and temperature-dependent `Tau_y`
- Yield strength estimates:
  - `YS 25-273C PRIOR`, `YS 25C PRIOR`, `YS T C PRIOR`, `YS T-273C PRIOR`
- Vickers hardness estimates from `HV = 3000 * Tau_y * (3/9.807) + 150`

## Outputs
- Results are written to `CalcFiles/STOIC_OUT_<index>.csv`.
- A tracker log is appended to `PROP-tracker.txt`.

## Gaps / Next Steps
- Provide `SPACE_v0.par` and `strength_model.py` to run the workflow.
- Validate the hard-coded property tables (elastic constants and elemental properties).
- Confirm test temperature strategy (currently fixed at 1300 C for all compositions).
- Confirm whether the element set should be expanded or constrained to ULTIMATE specs.

