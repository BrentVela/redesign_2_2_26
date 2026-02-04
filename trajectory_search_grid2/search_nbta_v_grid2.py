#!/usr/bin/env python
"""
2D grid search varying V and Ta, with Nb as balance, for:
- Cr=5 at%, W=10 at%, Zr=3 at%
- V <= 20 at%
- Target density 9.5 g/cc (ROM) and liquidus ~2750C (TC)

Outputs ranked candidates to CSV.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from tc_python import TCPython, CompositionUnit

ATOMIC_WEIGHTS = {
    "Nb": 92.906,
    "Ta": 180.948,
    "V": 50.942,
    "W": 183.84,
    "Cr": 51.996,
    "Zr": 91.224,
}

ELEMENT_DENSITIES = {
    "Nb": 8.57,
    "Ta": 16.69,
    "V": 6.11,
    "W": 19.25,
    "Cr": 7.15,
    "Zr": 6.52,
}


def density_rule_of_mixtures(comp):
    molar_mass = sum(comp[el] * ATOMIC_WEIGHTS[el] for el in comp)
    molar_volume = sum(comp[el] * ATOMIC_WEIGHTS[el] / ELEMENT_DENSITIES[el] for el in comp)
    return molar_mass / molar_volume


def solidus_tc_batch(comps):
    active = ["Nb", "Ta", "V", "Cr", "W", "Zr"]
    solidus = []
    with TCPython() as session:
        session.disable_caching()
        calc = (
            session
            .select_database_and_elements('TCHEA7', active)
            .get_system()
            .with_property_model_calculation("Liquidus and Solidus Temperature")
            .set_argument('upperTemperatureLimit', 4000)
            .set_composition_unit(CompositionUnit.MOLE_PERCENT)
            .set_temperature(300)
        )
        for comp in comps:
            for el in active[:-1]:
                calc.set_composition(el, comp[el] * 100.0)
            result = calc.calculate()
            solidus.append(result.get_value_of('Solidus temperature'))
    return solidus


def main():
    target_density = 9.5
    target_tm_c = 2750
    fixed = {"Cr": 0.05, "W": 0.10, "Zr": 0.03}

    # Grid over V and Ta; Nb is balance
    v_grid = np.linspace(0.0, 0.20, 21)
    ta_grid = np.linspace(0.0, 0.30, 31)

    rows = []
    for v in v_grid:
        for ta in ta_grid:
            remaining = 1.0 - fixed['Cr'] - fixed['W'] - fixed['Zr'] - v - ta
            if remaining < 0:
                continue
            nb = remaining
            comp = {"Cr": fixed['Cr'], "W": fixed['W'], "Zr": fixed['Zr'], "V": v, "Ta": ta, "Nb": nb}
            dens = density_rule_of_mixtures(comp)
            rows.append({**comp, 'Density_g_cc_ROM': dens})

    df = pd.DataFrame(rows)
    df['density_err'] = (df['Density_g_cc_ROM'] - target_density).abs()

    # Limit TC calculations to top candidates by density
    top = df.sort_values('density_err').head(120).reset_index(drop=True)

    comps = [row for row in top[['Nb','Ta','V','Cr','W','Zr']].to_dict('records')]
    tms = solidus_tc_batch(comps)
    top['Solidus_K_TC'] = tms
    top['Solidus_C_TC'] = top['Solidus_K_TC'] - 273.15
    top['tm_err'] = (top['Solidus_C_TC'] - target_tm_c).abs()

    # rank by density error then tm error
    top = top.sort_values(['density_err','tm_err']).reset_index(drop=True)

    out = Path('trajectory_search_grid2/search_results.csv')
    top.to_csv(out, index=False)
    print(f"Wrote {out}")

    best_density = top.iloc[0]
    best_tm = top.sort_values('tm_err').iloc[0]

    print('\nBest by density (ROM):')
    print(best_density[['Nb','Ta','V','Cr','W','Zr','Density_g_cc_ROM','Solidus_C_TC']].to_string())

    print('\nBest by liquidus (TC):')
    print(best_tm[['Nb','Ta','V','Cr','W','Zr','Density_g_cc_ROM','Solidus_C_TC']].to_string())


if __name__ == '__main__':
    main()
