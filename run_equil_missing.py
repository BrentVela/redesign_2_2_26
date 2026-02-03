import os
import numpy as np
import pandas as pd
from tc_python import TCPython, LoggingPolicy, ThermodynamicQuantity


def _ensure_tc_home() -> None:
    if os.environ.get("TC25B_HOME"):
        return
    default_path = "/opt/Thermo-Calc/2025b"
    if os.path.isdir(default_path):
        os.environ["TC25B_HOME"] = default_path
    else:
        raise EnvironmentError("Environment variable 'TC25B_HOME' not found and default path is missing.")


_ensure_tc_home()

MISSING_IDS = ["NbTaV65_path_04", "NbTaV65_path_09"]
TEMPS_C = [500, 600, 1300, 2000]
ELEMENTS = sorted(["Cr","Mo","Nb","Ta","Ti","V","W"])

results_df = pd.read_csv("trajectory.csv")
subset = results_df[results_df["Alloy_ID"].isin(MISSING_IDS)].copy()
if subset.empty:
    raise SystemExit("No matching alloys found")

rows = []
with TCPython(logging_policy=LoggingPolicy.NONE) as session:
    session.disable_caching()
    for _, row in subset.iterrows():
        comp = np.array([row[e] for e in ELEMENTS])
        active_el = [el for el, v in zip(ELEMENTS, comp) if v > 0]
        system = session.select_database_and_elements('TCHEA7', active_el).get_system()
        eq = system.with_single_equilibrium_calculation().set_condition(
            ThermodynamicQuantity.temperature(), 298
        )
        if len(active_el) == 1:
            eq.set_dependent_element(active_el[0])
        else:
            for j in range(len(active_el) - 1):
                eq.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(active_el[j]), comp[j])

        out = {
            "Alloy_ID": row["Alloy_ID"],
            **{el: row[el] for el in ELEMENTS},
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
