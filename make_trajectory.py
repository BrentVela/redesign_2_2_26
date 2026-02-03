#!/usr/bin/env python
import math
import pandas as pd

ATOMIC_WEIGHTS = {
    "Nb": 92.906,
    "Ta": 180.948,
    "V": 50.942,
    "Mo": 95.95,
    "W": 183.84,
    "Ti": 47.867,
    "Cr": 51.996,
}

ELEMENT_DENSITIES = {
    "Nb": 8.57,
    "Ta": 16.69,
    "V": 6.11,
    "Mo": 10.28,
    "W": 19.25,
    "Ti": 4.51,
    "Cr": 7.15,
}


def alloy_density(comp):
    molar_mass = sum(comp[el] * ATOMIC_WEIGHTS[el] for el in comp)
    molar_volume = sum(comp[el] * ATOMIC_WEIGHTS[el] / ELEMENT_DENSITIES[el] for el in comp)
    return molar_mass / molar_volume


def solve_v_for_density(target_density, tol=1e-9):
    remaining = 0.65
    lo, hi = 0.0, remaining
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        x_nb = (remaining - mid) / 2.0
        x_ta = x_nb
        comp = {
            "Nb": x_nb,
            "Ta": x_ta,
            "V": mid,
            "Mo": 0.10,
            "W": 0.10,
            "Ti": 0.10,
            "Cr": 0.05,
        }
        rho = alloy_density(comp)
        if rho > target_density:
            lo = mid
        else:
            hi = mid
        if abs(rho - target_density) < tol:
            break
    x_v = 0.5 * (lo + hi)
    return x_v


if __name__ == "__main__":
    # Baseline (Nb=Ta=V) for (NbTaV)65
    v_start = 0.65 / 3.0
    v_target = solve_v_for_density(9.5)

    steps = 10
    rows = []
    for i in range(steps + 1):
        frac = i / steps
        v = v_start + (v_target - v_start) * frac
        nb = (0.65 - v) / 2.0
        ta = nb
        comp = {
            "Cr": 0.05,
            "Mo": 0.10,
            "Nb": nb,
            "Ta": ta,
            "Ti": 0.10,
            "V": v,
            "W": 0.10,
        }
        density = alloy_density(comp)
        row = {"Alloy_ID": f"NbTaV65_path_{i:02d}", **comp, "Density_g_cc": density}
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("trajectory.csv", index=False)
    df.drop(columns=["Alloy_ID", "Density_g_cc"]).to_parquet("SPACE_v0.par", index=False)

    print("Wrote trajectory.csv")
    print("Wrote SPACE_v0.par")
    print(f"V_start={v_start:.12f} V_target={v_target:.12f}")
