import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Times", "Computer Modern Roman"],
    "font.size": 10,
    "axes.labelsize": 10,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "axes.linewidth": 1.0,
    "axes.formatter.useoffset": False,
    "axes.xmargin": 0,
    "axes.ymargin": 0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.format": "pdf",
    "savefig.transparent": False,
    "figure.constrained_layout.use": True,
    "figure.subplot.left": 0.125,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.2,
    "figure.subplot.top": 0.95,
    "axes.autolimit_mode": "round_numbers",
}


def apply_affine_style() -> None:
    mpl.rcParams.update(_PARAMS)


def main() -> None:
    df = pd.read_csv("CalcFiles/STOIC_EQUIL_OUT.csv")
    if df.empty:
        raise SystemExit("No data found in CalcFiles/STOIC_EQUIL_OUT.csv")

    bcc_cols = [
        "EQ 500C BCC_B2 MOL",
        "EQ 500C BCC_B2#2 MOL",
    ]
    for col in bcc_cols:
        if col not in df.columns:
            raise SystemExit(f"Missing column: {col}")

    bcc_sum = df[bcc_cols].sum(axis=1, skipna=True)

    apply_affine_style()
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    ax.plot(df["V"], bcc_sum, marker="o")
    ax.set_xlabel("V (at. frac.)")
    ax.set_ylabel("BCC fraction (500C)")
    ax.set_ylim(0, 1.05)
    plt.savefig("bcc_fraction_vs_v.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
