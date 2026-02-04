#!/usr/bin/env python
"""
Search compositions for (Nb, Ta, V, Cr, W, Zr) with constraints:
- Start from (NbTaV)65Cr5Mo10Ti10W10 baseline
- Remove Mo and Ti, keep Cr=5, W=10
- Add Zr=3 at%
- V <= 20 at%
- Adjust Ta to hit density 9.5 g/cc and liquidus ~2750 C (TC property model)
- If both can't be matched, prioritize density by adjusting Nb; maximize liquidus.

Outputs a CSV with candidates and the best-by-density and best-by-liquidus rankings.
"""

import math
from pathlib import Path
import pandas as pd
import numpy as np
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


def liquidus_solidus_tc_batch(comps):
    active = ["Nb", "Ta", "V", "Cr", "W", "Zr"]
    liquidus = []
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
            liquidus.append(result.get_value_of('Liquidus temperature'))
            solidus.append(result.get_value_of('Solidus temperature'))
    return liquidus, solidus


def main():
    target_density = 9.5
    target_tm_c = 2750

    fixed = {"Cr": 0.05, "W": 0.10, "Zr": 0.03}

    # sweep V from 0 to 0.20 (coarse)
    v_grid = np.linspace(0.0, 0.20, 21)

    rows = []
    for v in v_grid:
        remaining = 1.0 - fixed['Cr'] - fixed['W'] - fixed['Zr'] - v
        if remaining <= 0:
            continue

        # Solve Ta for density target using ROM bisection; Nb is balance.
        def dens_for_ta(ta):
            nb = remaining - ta
            comp = {"Cr": fixed['Cr'], "W": fixed['W'], "Zr": fixed['Zr'], "V": v, "Ta": ta, "Nb": nb}
            return density_rule_of_mixtures(comp)

        lo, hi = 0.0, remaining
        dens_lo = dens_for_ta(lo)
        dens_hi = dens_for_ta(hi)

        if (dens_lo - target_density) * (dens_hi - target_density) <= 0:
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                dens_mid = dens_for_ta(mid)
                if (dens_lo - target_density) * (dens_mid - target_density) <= 0:
                    hi = mid
                    dens_hi = dens_mid
                else:
                    lo = mid
                    dens_lo = dens_mid
            ta = 0.5 * (lo + hi)
            dens = dens_for_ta(ta)
        else:
            if abs(dens_lo - target_density) <= abs(dens_hi - target_density):
                ta = lo
                dens = dens_lo
            else:
                ta = hi
                dens = dens_hi

        nb = remaining - ta
        rows.append({
            "Nb": nb, "Ta": ta, "V": v, "Cr": fixed['Cr'], "W": fixed['W'], "Zr": fixed['Zr'],
            "Density_g_cc_ROM": dens,
        })

    df = pd.DataFrame(rows)
    df['density_err'] = (df['Density_g_cc_ROM'] - target_density).abs()

    # Compute TC liquidus for these candidates (reuse one TC session)
    comps = [row for row in df[['Nb','Ta','V','Cr','W','Zr']].to_dict('records')]
    liq, sol = liquidus_solidus_tc_batch(comps)
    df['Liquidus_K_TC'] = liq
    df['Solidus_K_TC'] = sol
    df['Liquidus_C_TC'] = df['Liquidus_K_TC'] - 273.15
    df['Solidus_C_TC'] = df['Solidus_K_TC'] - 273.15
    df['tm_err'] = (df['Liquidus_C_TC'] - target_tm_c).abs()
    df['ts_err'] = (df['Solidus_C_TC'] - target_tm_c).abs()

    # rank by density error then solidus error
    df = df.sort_values(['density_err','ts_err']).reset_index(drop=True)

    out = Path('trajectory_search/search_results.csv')
    df.to_csv(out, index=False)
    print(f"Wrote {out}")

    # summarize best by density and best by liquidus
    best_density = df.iloc[0]
    best_ts = df.sort_values('ts_err').iloc[0]

    print('\nBest by density (ROM):')
    print(best_density[['Nb','Ta','V','Cr','W','Zr','Density_g_cc_ROM','Liquidus_C_TC']].to_string())

    print('\nBest by solidus (TC):')
    print(best_ts[['Nb','Ta','V','Cr','W','Zr','Density_g_cc_ROM','Solidus_C_TC']].to_string())


if __name__ == '__main__':
    main()
