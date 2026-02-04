#!/usr/bin/env python
"""
Append RT (25C) THCD to search_results_solidus_density_with_props.csv.
Only runs THCD calculation to avoid re-running other properties.
"""
import os
import numpy as np
import pandas as pd
from tc_python import TCPython, LoggingPolicy, ThermodynamicQuantity


def ensure_tc_home() -> None:
    if os.environ.get("TC25B_HOME"):
        return
    default_path = "/opt/Thermo-Calc/2025b"
    if os.path.isdir(default_path):
        os.environ["TC25B_HOME"] = default_path
    else:
        raise EnvironmentError("TC25B_HOME not set and default path missing")


def main() -> None:
    ensure_tc_home()
    in_path = "trajectory_search/search_results_solidus_density_with_props.csv"
    df = pd.read_csv(in_path)
    elements = ["Cr", "Nb", "Ta", "V", "W", "Zr"]

    rt_thcd = []
    with TCPython(logging_policy=LoggingPolicy.NONE) as session:
        session.disable_caching()
        system = (
            session
            .select_database_and_elements('TCHEA7', elements)
            .get_system()
        )
        eq = system.with_single_equilibrium_calculation().set_condition(
            ThermodynamicQuantity.temperature(), 25 + 273.15
        )

        for _, row in df.iterrows():
            comp = np.array([row[e] for e in elements])
            for j in range(len(elements) - 1):
                eq.set_condition(
                    ThermodynamicQuantity.mole_fraction_of_a_component(elements[j]),
                    comp[j],
                )
            res = eq.calculate()
            try:
                rt_thcd.append(res.get_value_of('THCD'))
            except Exception:
                rt_thcd.append(np.nan)

    df['EQ 25C THCD (W/mK)'] = rt_thcd
    df.to_csv(in_path, index=False)
    print(f"Wrote {in_path} with EQ 25C THCD (W/mK)")


if __name__ == "__main__":
    main()
