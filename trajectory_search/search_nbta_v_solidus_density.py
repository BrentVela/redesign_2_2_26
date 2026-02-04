#!/usr/bin/env python
"""
Search compositions for (Nb, Ta, V, Cr, W, Zr) with constraints:
- Remove Mo and Ti, keep Cr=5, W=10
- Add Zr=3 at%
- V <= 20 at% (1 at% step)
- Target solidus 2750C; then minimize density

Method: for each V, solve Ta by bisection to match solidus.
Outputs CSV with 21 candidates.
"""

import numpy as np
import pandas as pd
from pathlib import Path
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


def main():
    target_solidus_c = 2750

    fixed = {"Cr": 0.05, "W": 0.10, "Zr": 0.03}
    v_grid = np.linspace(0.0, 0.20, 21)

    rows = []
    with TCPython() as session:
        session.disable_caching()
        calc = (
            session
            .select_database_and_elements('TCHEA7', ["Nb", "Ta", "V", "Cr", "W", "Zr"])
            .get_system()
            .with_property_model_calculation("Liquidus and Solidus Temperature")
            .set_argument('upperTemperatureLimit', 4000)
            .set_composition_unit(CompositionUnit.MOLE_PERCENT)
            .set_temperature(300)
        )

        for v in v_grid:
            remaining = 1.0 - fixed['Cr'] - fixed['W'] - fixed['Zr'] - v
            if remaining <= 0:
                continue

            def eval_ta(ta):
                nb = remaining - ta
                comp = {"Cr": fixed['Cr'], "W": fixed['W'], "Zr": fixed['Zr'], "V": v, "Ta": ta, "Nb": nb}
                for el in ["Nb", "Ta", "V", "Cr", "W"]:
                    calc.set_composition(el, comp[el] * 100.0)
                result = calc.calculate()
                solidus_c = result.get_value_of('Solidus temperature') - 273.15
                liquidus_c = result.get_value_of('Liquidus temperature') - 273.15
                dens = density_rule_of_mixtures(comp)
                return solidus_c, liquidus_c, dens

            lo, hi = 0.0, remaining
            sol_lo, liq_lo, dens_lo = eval_ta(lo)
            sol_hi, liq_hi, dens_hi = eval_ta(hi)

            # If target is bracketed, bisection; otherwise choose closest end
            if (sol_lo - target_solidus_c) * (sol_hi - target_solidus_c) <= 0:
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    sol_mid, liq_mid, dens_mid = eval_ta(mid)
                    if (sol_lo - target_solidus_c) * (sol_mid - target_solidus_c) <= 0:
                        hi, sol_hi, liq_hi, dens_hi = mid, sol_mid, liq_mid, dens_mid
                    else:
                        lo, sol_lo, liq_lo, dens_lo = mid, sol_mid, liq_mid, dens_mid
                ta = 0.5 * (lo + hi)
                sol_c, liq_c, dens = eval_ta(ta)
            else:
                if abs(sol_lo - target_solidus_c) <= abs(sol_hi - target_solidus_c):
                    ta, sol_c, liq_c, dens = lo, sol_lo, liq_lo, dens_lo
                else:
                    ta, sol_c, liq_c, dens = hi, sol_hi, liq_hi, dens_hi

            nb = remaining - ta
            rows.append({
                "Nb": nb, "Ta": ta, "V": v, "Cr": fixed['Cr'], "W": fixed['W'], "Zr": fixed['Zr'],
                "Density_g_cc_ROM": dens,
                "Solidus_C_TC": sol_c,
                "Liquidus_C_TC": liq_c,
            })

    df = pd.DataFrame(rows)
    df['solidus_err'] = (df['Solidus_C_TC'] - target_solidus_c).abs()

    # Rank by solidus error, then lowest density
    df = df.sort_values(['solidus_err','Density_g_cc_ROM']).reset_index(drop=True)

    out = Path('trajectory_search/search_results_solidus_density.csv')
    df.to_csv(out, index=False)
    print(f"Wrote {out}")

    best = df.iloc[0]
    print('\nBest by solidus then density:')
    print(best[['Nb','Ta','V','Cr','W','Zr','Density_g_cc_ROM','Solidus_C_TC','Liquidus_C_TC']].to_string())


if __name__ == '__main__':
    main()
