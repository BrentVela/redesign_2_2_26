import argparse
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_PARAMS = {
    # ---------- Fonts ----------
    "font.family": "serif",
    "font.serif": ["Times", "Computer Modern Roman"],
    "font.size": 10,
    "axes.labelsize": 10,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{bm} \usepackage{siunitx} \usepackage{chemformula}",

    # ---------- Lines & markers ----------
    "lines.linewidth": 2.0,
    "lines.markersize": 6,

    # ---------- Axes ----------
    "axes.linewidth": 1.0,
    "axes.labelweight": "normal",
    "axes.formatter.useoffset": False,
    "axes.xmargin": 0,
    "axes.ymargin": 0,

    # ---------- Ticks ----------
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

    # ---------- Figure & save ----------
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.format": "pdf",
    "savefig.transparent": False,

    # ---------- Legend ----------
    "legend.fontsize": 10,
    "legend.frameon": False,
    "legend.numpoints": 1,
    "legend.handlelength": 2,
    "legend.scatterpoints": 1,
    "legend.labelspacing": 0.5,
    "legend.markerscale": 0.9,
    "legend.handletextpad": 0.5,
    "legend.borderaxespad": 0.5,
    "legend.borderpad": 0.5,
    "legend.columnspacing": 1,

    # ---------- Subplot spacing ----------
    "figure.constrained_layout.use": True,
    "figure.subplot.left": 0.125,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.2,
    "figure.subplot.top": 0.95,
    "axes.autolimit_mode": "round_numbers",
}


_LABELS = {
    "V": r"V (at. frac.)",
    "index": r"Trajectory index",
    "Density Avg": r"Density (\si{\gram\per\centi\meter\cubed})",
    "Tm Avg": r"$T_m$ (\si{\kelvin})",
    "VEC Avg": r"VEC",
    "Bulk Modulus Avg": r"Bulk modulus (\si{\giga\pascal})",
    "Shear Modulus Avg": r"Shear modulus (\si{\giga\pascal})",
    "Cauchy Pres Avg": r"Cauchy pressure (\si{\giga\pascal})",
    "Pugh_Ratio_PRIOR": r"Pugh ratio",
    "Zener Ratio": r"Zener ratio",
    "Universal Anisotropy": r"Universal anisotropy",
    "YS 25C PRIOR": r"Yield strength (\si{\mega\pascal})",
    "YS T C PRIOR": r"Yield strength (\si{\mega\pascal})",
    "HV 25C PRIOR": r"Vickers hardness",
    "HV T C PRIOR": r"Vickers hardness",
}


def apply_affine_style(use_tex: bool = True) -> None:
    params = dict(_PARAMS)
    if not use_tex:
        params["text.usetex"] = False
    mpl.rcParams.update(params)


def _label_for(col: str) -> str:
    return _LABELS.get(col, col)


def _parse_props(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def plot_properties(
    df: pd.DataFrame,
    xcol: str,
    props: list[str],
    out_path: str,
    use_tex: bool,
) -> None:
    apply_affine_style(use_tex=use_tex)

    n = len(props)
    if n == 0:
        raise ValueError("No properties specified for plotting.")

    ncols = 2 if n > 1 else 1
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6.8, 2.6 * nrows))

    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = [axes]

    x = df[xcol]
    for ax, prop in zip(axes, props):
        if prop not in df.columns:
            raise ValueError(f"Property '{prop}' not found in CSV.")
        ax.plot(x, df[prop], marker="o")
        ax.set_xlabel(_label_for(xcol))
        ax.set_ylabel(_label_for(prop))

    for ax in axes[len(props):]:
        ax.axis("off")

    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot properties vs trajectory using affine projection styling.")
    parser.add_argument("--csv", default="CalcFiles/STOIC_OUT_0.csv")
    parser.add_argument("--xcol", default="V")
    parser.add_argument(
        "--props",
        default=(
            "Density Avg,Tm Avg,VEC Avg,Cauchy Pres Avg,Pugh_Ratio_PRIOR,YS 25C PRIOR"
        ),
    )
    parser.add_argument("--out", default="properties_vs_trajectory.pdf")
    parser.add_argument("--no-tex", action="store_true", help="Disable LaTeX rendering.")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if df.empty:
        raise ValueError("No data found in the CSV.")
    if args.xcol not in df.columns:
        raise ValueError(f"x-column '{args.xcol}' not found in CSV.")

    props = _parse_props(args.props)
    plot_properties(df, args.xcol, props, args.out, use_tex=not args.no_tex)


if __name__ == "__main__":
    main()
