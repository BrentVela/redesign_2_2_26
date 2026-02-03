import pandas as pd
import numpy as np
from tc_python import TCPython, LoggingPolicy, ThermodynamicQuantity, SingleEquilibriumOptions, ALL_PHASES

TEMPS_C = [500, 600, 1300, 2000]
ELEMENTS = sorted(["Cr", "Mo", "Nb", "Ta", "Ti", "V", "W"])

# Load existing EQUIL outputs
files = ["CalcFiles/EQUIL_AP_OUT_0.csv", "CalcFiles/EQUIL_AP_OUT_5.csv"]

equil_frames = []
for path in files:
    df = pd.read_csv(path)
    unnamed = [c for c in df.columns if c.lower().startswith("unnamed")]
    if unnamed:
        df = df.drop(columns=unnamed)
    df["__source_file"] = path
    equil_frames.append(df)

equil = pd.concat(equil_frames, ignore_index=True)

if "EQ 500C BCC_B2 MOL" not in equil.columns:
    raise SystemExit("Missing EQ 500C BCC_B2 MOL column")

missing = equil[equil["EQ 500C BCC_B2 MOL"].isna()].copy()
missing_ids = missing["Alloy_ID"].dropna().unique().tolist()
if not missing_ids:
    print("No missing alloys found")
    raise SystemExit(0)

traj = pd.read_csv("trajectory.csv")
subset = traj[traj["Alloy_ID"].isin(missing_ids)].copy()
if subset.empty:
    raise SystemExit("No matching alloys in trajectory.csv")

rows = []
with TCPython(logging_policy=LoggingPolicy.NONE) as session:
    session.disable_caching()
    options = (
        SingleEquilibriumOptions()
        .set_required_accuracy(1.0e-2)
        .set_max_no_of_iterations(300)
        .set_smallest_fraction(1.0e-6)
    )
    for _, row in subset.iterrows():
        comp = np.array([row[e] for e in ELEMENTS])
        active_el = [el for el, v in zip(ELEMENTS, comp) if v > 0]
        system = session.select_database_and_elements('TCHEA7', active_el).get_system()
    eq = (
        system.with_single_equilibrium_calculation()
        .set_condition(ThermodynamicQuantity.temperature(), 298)
        .with_options(options)
    )
    eq.disable_global_minimization()
    eq.set_phase_to_suspended(ALL_PHASES)
    for phase in ["BCC_B2", "BCC_B2#2", "C14_LAVES", "C15_LAVES", "C36_LAVES"]:
        eq.set_phase_to_entered(phase)
        if len(active_el) == 1:
            eq.set_dependent_element(active_el[0])
        else:
            for j in range(len(active_el) - 1):
                eq.set_condition(
                    ThermodynamicQuantity.mole_fraction_of_a_component(active_el[j]),
                    comp[j],
                )

        out = {
            "Alloy_ID": row["Alloy_ID"],
            **{el: row[el] for el in ELEMENTS},
            "Density_g_cc": row.get("Density_g_cc", np.nan),
        }
        for temp in TEMPS_C:
            try:
                eq_temp = eq.set_condition(ThermodynamicQuantity.temperature(), temp + 273)
                res = eq_temp.calculate()
                for phase in res.get_phases():
                    frac = res.get_value_of(f"NPM({phase})")
                    if frac > 0:
                        out[f"EQ {temp}C {phase} MOL"] = frac
            except Exception:
                # Leave missing entries as NaN
                pass
        rows.append(out)

out_df = pd.DataFrame(rows)
out_df.to_csv("CalcFiles/EQUIL_AP_OUT_MISSING.csv", index=False)
print("Wrote CalcFiles/EQUIL_AP_OUT_MISSING.csv")

# Merge new results into original EQUIL files
merged = equil.drop(columns=["__source_file"]).copy()
for _, new_row in out_df.iterrows():
    alloy_id = new_row["Alloy_ID"]
    for col, val in new_row.items():
        if col == "Alloy_ID":
            continue
        if col not in merged.columns:
            merged[col] = np.nan
        merged.loc[merged["Alloy_ID"] == alloy_id, col] = val

# Write back to original files (preserving row membership)
for path in files:
    part = merged[merged["Alloy_ID"].isin(equil[equil["__source_file"] == path]["Alloy_ID"])].copy()
    part.to_csv(path, index=False)

print("Updated EQUIL_AP_OUT files")
