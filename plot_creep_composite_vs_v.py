#!/usr/bin/env python
"""
Plot a priori composite creep rates (min of NH, Coble, and PL5) vs V content.

Reads precomputed rates from CalcFiles/CREEP_OUT_0.csv and compositions from
trajectory.csv. Outputs a CSV and a PDF plot for a chosen temperature.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--temp', type=int, default=1200, help='Temperature in C')
    ap.add_argument('--out', default=None, help='Output PDF path')
    args = ap.parse_args()

    creep_path = Path('CalcFiles/CREEP_OUT_0.csv')
    traj_path = Path('trajectory.csv')
    if not creep_path.exists():
        raise SystemExit(f"Missing {creep_path}")
    if not traj_path.exists():
        raise SystemExit(f"Missing {traj_path}")

    creep = pd.read_csv(creep_path)
    traj = pd.read_csv(traj_path)

    t = args.temp
    col_nh = f"{t} Min Creep NH [1/s]"
    col_cb = f"{t} Min Creep CB [1/s]"
    col_pl = f"{t} Min Creep PL5 [1/s]"

    for col in (col_nh, col_cb, col_pl):
        if col not in creep.columns:
            raise SystemExit(f"Missing column: {col}")

    out = traj.copy()
    out[col_nh] = creep[col_nh]
    out[col_cb] = creep[col_cb]
    out[col_pl] = creep[col_pl]

    # Composite a priori rate: minimum of NH, Coble, PL5
    out[f"{t} Min Creep Composite [1/s]"] = out[[col_nh, col_cb, col_pl]].min(axis=1)

    # Save CSV
    csv_out = Path(f'CalcFiles/creep_composite_{t}C_vs_v.csv')
    out.to_csv(csv_out, index=False)

    # Plot vs V
    if 'V' not in out.columns:
        raise SystemExit("trajectory.csv missing V column")

    x = out['V']
    y = out[f"{t} Min Creep Composite [1/s]"]

    fig, ax = plt.subplots(figsize=(6,4), dpi=200)
    ax.scatter(x, np.log10(y), s=30, color='#1a9850', edgecolor='black', linewidth=0.4)
    ax.set_xlabel('V (at. fraction)')
    ax.set_ylabel('Composite creep rate, log10(1/s)')
    ax.set_title(f'Composite creep rate vs V at {t}C')
    ax.grid(True, linestyle=':', linewidth=0.6, alpha=0.6)
    fig.tight_layout()

    pdf_out = Path(args.out) if args.out else Path(f'creep_composite_{t}C_vs_v.pdf')
    fig.savefig(pdf_out)

    print(f"Wrote {csv_out}")
    print(f"Wrote {pdf_out}")


if __name__ == '__main__':
    main()
