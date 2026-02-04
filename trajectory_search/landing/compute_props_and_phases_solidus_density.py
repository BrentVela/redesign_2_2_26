#!/usr/bin/env python
"""
Compute derived properties (VEC, Pugh ratio, Cauchy pressure, YS@1300C) and
Thermo-Calc phase fractions + thermal conductivity at 600C and 1300C
for candidates in trajectory_search/landing/search_results_solidus_density.csv.
"""
import os
import sys
import numpy as np
import pandas as pd
from tc_python import TCPython, LoggingPolicy, ThermodynamicQuantity
from tc_python.exceptions import CalculationException

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import strength_model


def ensure_tc_home() -> None:
    if os.environ.get("TC25B_HOME"):
        return
    default_path = "/opt/Thermo-Calc/2025b"
    if os.path.isdir(default_path):
        os.environ["TC25B_HOME"] = default_path
    else:
        raise EnvironmentError("TC25B_HOME not set and default path missing")


# Elemental data (consistent with 1_find_stoich_props_parallel_for_Class.py)
PROP_DATA = {
    'Cr': {'Valence Electrons': 6},
    'Mo': {'Valence Electrons': 6},
    'Nb': {'Valence Electrons': 5},
    'Ta': {'Valence Electrons': 5},
    'Ti': {'Valence Electrons': 4},
    'V':  {'Valence Electrons': 5},
    'W':  {'Valence Electrons': 6},
    'Zr': {'Valence Electrons': 4},
}

ELAST_DATA = {
    'Cr': {'C11': 247.6, 'C12': 73.4,  'C44': 48.3,  'B': 131.5, 'G': 48.3},
    'Mo': {'C11': 466.0, 'C12': 165.2, 'C44': 99.5,  'B': 265.8, 'G': 99.5},
    'Nb': {'C11': 247.2, 'C12': 140.0, 'C44': 14.2,  'B': 175.7, 'G': 14.2},
    'Ta': {'C11': 260.9, 'C12': 165.2, 'C44': 70.4,  'B': 197.1, 'G': 70.4},
    'Ti': {'C11': 95.9,  'C12': 115.9, 'C44': 40.3,  'B': 109.2, 'G': 40.3},
    'V':  {'C11': 272.0, 'C12': 144.8, 'C44': 17.6,  'B': 187.2, 'G': 17.6},
    'W':  {'C11': 517.8, 'C12': 201.7, 'C44': 139.4, 'B': 306.4, 'G': 139.4},
    'Zr': {'C11': 81.8,  'C12': 94.3,  'C44': 30.2,  'B': 90.1,  'G': 30.2},
}


def compute_basic_props(df: pd.DataFrame, elements: list[str]) -> pd.DataFrame:
    comp = df[elements]
    el_c11 = np.array([ELAST_DATA[e]['C11'] for e in elements])
    el_c12 = np.array([ELAST_DATA[e]['C12'] for e in elements])
    el_c44 = np.array([ELAST_DATA[e]['C44'] for e in elements])
    el_b = np.array([ELAST_DATA[e]['B'] for e in elements])
    el_g = np.array([ELAST_DATA[e]['G'] for e in elements])
    el_vec = np.array([PROP_DATA[e]['Valence Electrons'] for e in elements])

    c11 = comp.values @ el_c11
    c12 = comp.values @ el_c12
    c44 = comp.values @ el_c44
    b_avg = comp.values @ el_b
    g_avg = comp.values @ el_g
    vec = comp.values @ el_vec

    cauchy = c12 - c44
    pugh = b_avg / g_avg

    ys_1300 = []
    for _, row in df.iterrows():
        comp_i = np.array([row[e] for e in elements])
        _, _, tau = strength_model.model_Control(elements, comp_i, T=1300 + 273.15)
        ys_1300.append(3000 * tau)

    df = df.copy()
    df['VEC Avg'] = vec
    df['Cauchy Pres Avg'] = cauchy
    df['Pugh_Ratio_PRIOR'] = pugh
    df['YS_1300C_MPa'] = ys_1300
    return df


def tc_phases_and_thcd(df: pd.DataFrame, elements: list[str]) -> pd.DataFrame:
    temps_c = [600, 1300]
    rows = []
    with TCPython(logging_policy=LoggingPolicy.NONE) as session:
        session.disable_caching()
        system = (
            session
            .select_database_and_elements('TCHEA7', elements)
            .get_system()
        )
        eq = system.with_single_equilibrium_calculation().set_condition(
            ThermodynamicQuantity.temperature(), 298
        )

        for _, row in df.iterrows():
            comp = np.array([row[e] for e in elements])
            for j in range(len(elements) - 1):
                eq.set_condition(
                    ThermodynamicQuantity.mole_fraction_of_a_component(elements[j]),
                    comp[j],
                )

            out = {"Alloy_ID": row.get("Alloy_ID", None)}
            for temp in temps_c:
                eq_t = eq.set_condition(ThermodynamicQuantity.temperature(), temp + 273.15)
                try:
                    res = eq_t.calculate()
                except CalculationException:
                    continue
                else:
                    for phase in res.get_phases():
                        frac = res.get_value_of(f"NPM({phase})")
                        if frac > 0:
                            out[f"EQ {temp}C {phase} MOL"] = frac
                    try:
                        out[f"EQ {temp}C THCD (W/mK)"] = res.get_value_of('THCD')
                    except Exception:
                        pass

            rows.append(out)

    return pd.DataFrame(rows)


def main() -> None:
    ensure_tc_home()
    in_path = "trajectory_search/landing/search_results_solidus_density.csv"
    out_path = "trajectory_search/landing/search_results_solidus_density_with_props.csv"

    df = pd.read_csv(in_path).reset_index(drop=True)
    df["Alloy_ID"] = [f"TSD_{i:03d}" for i in range(len(df))]
    elements = ["Cr", "Nb", "Ta", "V", "W", "Zr"]

    df = compute_basic_props(df, elements)
    tc_df = tc_phases_and_thcd(df, elements)

    merged = pd.merge(df, tc_df, on="Alloy_ID", how="left")

    for temp in [600, 1300]:
        bcc_cols = [c for c in [f"EQ {temp}C BCC_B2 MOL", f"EQ {temp}C BCC_B2#2 MOL"] if c in merged.columns]
        laves_cols = [c for c in [f"EQ {temp}C C14_LAVES MOL", f"EQ {temp}C C15_LAVES MOL", f"EQ {temp}C C36_LAVES MOL"] if c in merged.columns]
        if bcc_cols:
            merged[f"EQ {temp}C BCC_SUM MOL"] = merged[bcc_cols].sum(axis=1, skipna=True)
        if laves_cols:
            merged[f"EQ {temp}C LAVES_SUM MOL"] = merged[laves_cols].sum(axis=1, skipna=True)

    merged.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
