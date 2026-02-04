import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import strength_model

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


def apply_affine_style() -> None:
    mpl.rcParams.update(_PARAMS)


def _format_comp_value(value: float) -> str:
    if np.isclose(value, round(value)):
        return str(int(round(value)))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _build_composition_label(row: pd.Series, elements: list[str]) -> str:
    parts = []
    for element in elements:
        value = row.get(element)
        if value is None or pd.isna(value) or value == 0:
            continue
        parts.append(f"{element}$_{{{_format_comp_value(float(value))}}}$")
    return "".join(parts)


def load_composition_map(csv_path: str) -> dict[int, str]:
    comp_df = pd.read_csv(csv_path)
    if "Alloys ID" not in comp_df.columns:
        return {}

    exclude = {
        "Iteration",
        "Alloys ID",
        "YS RT PRIOR (MPa)",
        "YS 600C PRIOR (MPa)",
        "YS 650C PRIOR (MPa)",
        "RT Compression YS (MPa)",
    }
    element_cols = [col for col in comp_df.columns if col not in exclude]

    comp_df = comp_df.dropna(subset=["Alloys ID"]).copy()
    comp_df["Alloys ID"] = pd.to_numeric(comp_df["Alloys ID"], errors="coerce").astype("Int64")

    comp_map: dict[int, str] = {}
    for _, row in comp_df.iterrows():
        alloy_id = row["Alloys ID"]
        if pd.isna(alloy_id):
            continue
        label = _build_composition_label(row, element_cols)
        if label:
            comp_map[int(alloy_id)] = label
    return comp_map


def load_compositions(csv_path: str) -> tuple[dict[int, list[float]], list[str]]:
    comp_df = pd.read_csv(csv_path)
    if "Alloys ID" not in comp_df.columns:
        return {}, []

    exclude = {
        "Iteration",
        "Alloys ID",
        "YS RT PRIOR (MPa)",
        "YS 600C PRIOR (MPa)",
        "YS 650C PRIOR (MPa)",
        "RT Compression YS (MPa)",
    }
    element_cols = [col for col in comp_df.columns if col not in exclude]

    comp_df = comp_df.dropna(subset=["Alloys ID"]).copy()
    comp_df["Alloys ID"] = pd.to_numeric(comp_df["Alloys ID"], errors="coerce").astype("Int64")

    compositions: dict[int, list[float]] = {}
    for _, row in comp_df.iterrows():
        alloy_id = row["Alloys ID"]
        if pd.isna(alloy_id):
            continue
        values = []
        for element in element_cols:
            value = row.get(element)
            if value is None or pd.isna(value):
                value = 0.0
            values.append(float(value) / 100.0)
        compositions[int(alloy_id)] = values
    return compositions, element_cols


def build_prediction_df(
    alloy_ids: list[int],
    compositions: dict[int, list[float]],
    elements: list[str],
    temp_min: float,
    temp_max: float,
    temp_step: float,
) -> pd.DataFrame:
    temps = np.arange(temp_min, temp_max + 0.0001, temp_step)
    rows = []
    for alloy_id in alloy_ids:
        comp = compositions.get(alloy_id)
        if comp is None:
            continue
        for temp_c in temps:
            _, _, tau = strength_model.model_Control(elements, comp, T=temp_c + 273.15)
            ys_mpa = 3000.0 * float(tau)
            rows.append({
                "Alloys ID": alloy_id,
                "Type": "Prediction",
                "Temp C": float(temp_c),
                "YS (MPa)": ys_mpa,
            })
    return pd.DataFrame(rows)


def _render_plot(df: pd.DataFrame, composition_map: dict[int, str]) -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    ax.set_title("Compressive Yield Strength Experiment vs Predictions")

    alloys = sorted(df["Alloys ID"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(alloys), 1)))
    color_map = {alloy: colors[i % len(colors)] for i, alloy in enumerate(alloys)}

    type_lower = df["Type"].fillna("").str.lower()
    for alloy in alloys:
        alloy_df = df[df["Alloys ID"] == alloy].sort_values("Temp C")
        alloy_type = type_lower.loc[alloy_df.index]
        pred = alloy_df[alloy_type == "prediction"]
        exp = alloy_df[alloy_type == "experiment"]

        if not pred.empty:
            ax.plot(pred["Temp C"], pred["YS (MPa)"], color=color_map[alloy], alpha=0.85)
        if not exp.empty:
            ax.scatter(
                exp["Temp C"],
                exp["YS (MPa)"],
                color=color_map[alloy],
                marker="o",
                s=30,
                edgecolors="none",
                alpha=0.9,
            )

    ax.set_xlabel(r"Temperature (\si{\celsius})")
    ax.set_ylabel(r"Yield strength (\si{\mega\pascal})")
    ax.set_xlim(0, 700)
    ax.axhline(250, color="gray", linestyle="--", linewidth=1.0, alpha=0.8)
    ax.axvline(650, color="gray", linestyle="--", linewidth=1.0, alpha=0.8)

    legend_items = []
    for alloy in alloys:
        label = composition_map.get(int(alloy), f"Alloy {int(alloy)}")
        legend_items.append(Line2D([0], [0], color=color_map[alloy], lw=2, label=label))
    if legend_items:
        ax.legend(handles=legend_items, loc="best", fontsize=7)

    ax.text(200, 1300, "MS5.2 Tensile Yield Strength", fontsize=7, color="black",
            ha="left", va="center")
    ax.plot([140, 190], [1300, 1300], color="gray", linestyle="--", linewidth=1.0, alpha=0.8)

    return fig, ax


def plot_ys_vs_temp(df: pd.DataFrame, out_path: str, composition_map: dict[int, str]) -> None:
    apply_affine_style()
    fig, _ = _render_plot(df, composition_map)
    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot YS vs temperature with affine projection styling.")
    parser.add_argument("--csv", default="YS_vs_T_data_points_IT0_1.csv")
    parser.add_argument("--out", default="YS_vs_T_plot_IT0_1.pdf")
    parser.add_argument("--comp-csv", default="IT0_1_YS_prior.csv")
    parser.add_argument("--use-model", action="store_true", help="Generate prediction curves with strength_model.")
    parser.add_argument("--temp-min", type=float, default=0.0)
    parser.add_argument("--temp-max", type=float, default=700.0)
    parser.add_argument("--temp-step", type=float, default=25.0)
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if df.empty:
        raise ValueError("No data found in the CSV.")

    composition_map = load_composition_map(args.comp_csv)

    if args.use_model:
        compositions, elements = load_compositions(args.comp_csv)
        exp_df = df[df["Type"].fillna("").str.lower() == "experiment"]
        alloy_ids = sorted(exp_df["Alloys ID"].dropna().unique()) if not exp_df.empty else sorted(compositions.keys())
        pred_df = build_prediction_df(alloy_ids, compositions, elements, args.temp_min, args.temp_max, args.temp_step)
        df = pd.concat([pred_df, exp_df], ignore_index=True)

    plot_ys_vs_temp(df, args.out, composition_map)


if __name__ == "__main__":
    main()
