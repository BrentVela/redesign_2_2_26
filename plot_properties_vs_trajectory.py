import argparse
import math
import os
import tempfile

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
    "text.usetex": False,

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
    "Density Avg": "Density (g/cc)",
    "Tm Avg": "Tₘ (K)",
    "VEC Avg": "VEC",
    "Bulk Modulus Avg": "Bulk modulus (GPa)",
    "Shear Modulus Avg": "Shear modulus (GPa)",
    "Cauchy Pres Avg": "Cauchy pressure (GPa)",
    "Pugh_Ratio_PRIOR": "Pugh ratio",
    "Zener Ratio": "Zener ratio",
    "Universal Anisotropy": "Universal anisotropy",
    "YS 25C PRIOR": "YS₂₅°C (MPa)",
    "YS T C PRIOR": "YSₜ°C (MPa)",
    "HV 25C PRIOR": "HV₂₅C",
    "HV T C PRIOR": "HVₜC",
}


def apply_affine_style(use_tex: bool = False) -> None:
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
    use_constrained: bool,
    single_column: bool,
) -> None:
    apply_affine_style(use_tex=use_tex)

    n = len(props)
    if n == 0:
        raise ValueError("No properties specified for plotting.")

    ncols = 1 if single_column else (2 if n > 1 else 1)
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(6.2 if ncols == 1 else 6.8, 2.6 * nrows),
        constrained_layout=use_constrained,
    )

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

    if not use_constrained:
        fig.tight_layout()

    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def prewarm_tex_cache(use_constrained: bool) -> None:
    apply_affine_style(use_tex=True)
    fig, ax = plt.subplots(figsize=(2.0, 1.5), constrained_layout=use_constrained)
    ax.plot([0, 1], [0, 1])
    ax.set_xlabel(r"Temperature (\si{\celsius})")
    ax.set_ylabel(r"Yield strength (\si{\mega\pascal})")
    if not use_constrained:
        fig.tight_layout()
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=True) as tmp:
        fig.savefig(tmp.name, bbox_inches="tight")
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
    parser.add_argument("--format", choices=["pdf", "png", "svg"], help="Override output format.")
    parser.add_argument("--no-constrained", action="store_true", help="Disable constrained_layout for speed.")
    parser.add_argument("--one-per-figure", action="store_true", help="Write one property per figure.")
    parser.add_argument("--out-dir", default=".", help="Output directory when using one-per-figure.")
    parser.add_argument("--single-column", action="store_true", help="Use one column layout for multi-plot output.")
    parser.add_argument("--prewarm-tex", action="store_true", help="Prewarm LaTeX cache before plotting.")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if df.empty:
        raise ValueError("No data found in the CSV.")
    if args.xcol not in df.columns:
        raise ValueError(f"x-column '{args.xcol}' not found in CSV.")

    props = _parse_props(args.props)
    use_tex = not args.no_tex
    use_constrained = not args.no_constrained

    if args.prewarm_tex and use_tex:
        prewarm_tex_cache(use_constrained)

    out_path = args.out
    if args.format:
        stem = out_path.rsplit(".", 1)[0]
        out_path = f"{stem}.{args.format}"

    if args.one_per_figure:
        os.makedirs(args.out_dir, exist_ok=True)
        for prop in props:
            prop_slug = prop.replace(" ", "_").replace("/", "_")
            stem = os.path.join(args.out_dir, f"properties_vs_trajectory_{prop_slug}")
            prop_out = f"{stem}.{out_path.rsplit('.', 1)[-1]}"
            plot_properties(
                df,
                args.xcol,
                [prop],
                prop_out,
                use_tex=use_tex,
                use_constrained=use_constrained,
                single_column=True,
            )
    else:
        plot_properties(
            df,
            args.xcol,
            props,
            out_path,
            use_tex=use_tex,
            use_constrained=use_constrained,
            single_column=args.single_column,
        )


if __name__ == "__main__":
    main()
