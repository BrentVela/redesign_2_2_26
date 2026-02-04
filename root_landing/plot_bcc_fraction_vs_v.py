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


def _load_bcc_vs_v(path: str, total_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        raise SystemExit(f"No data found in {path}")
    if "V" not in df.columns:
        raise SystemExit(f"Missing column 'V' in {path}")
    if total_col not in df.columns:
        raise SystemExit(f"Missing column '{total_col}' in {path}")
    return df[["V", total_col]]


def main() -> None:
    df_600 = _load_bcc_vs_v("CalcFiles/BCC_600C_fractions.csv", "BCC_600C_TOTAL_MOL")
    df_1300 = _load_bcc_vs_v("CalcFiles/BCC_1300C_fractions.csv", "BCC_1300C_TOTAL_MOL")

    apply_affine_style()
    fig, ax = plt.subplots(figsize=(4.4, 3.3))
    ax.plot(df_600["V"], df_600["BCC_600C_TOTAL_MOL"], marker="o", label="600C")
    ax.plot(df_1300["V"], df_1300["BCC_1300C_TOTAL_MOL"], marker="s", label="1300C")
    ax.set_xlabel("V (at. frac.)")
    ax.set_ylabel("BCC fraction")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False)
    plt.savefig("bcc_fraction_vs_v_600C_1300C.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
