#!/usr/bin/env python
"""
End-to-end featurization for Cafer's spreadsheets.
Adds: density (ROM), solidus/liquidus (TC), RT THCD (TC), YS@1300C,
Pugh ratio, Cauchy pressure, and VEC.
"""
import re
from datetime import datetime
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import importlib.util
import sys

from TC_Property_Module import compute_solidus_liquidus
from TC_Equilibrium_Module import compute_rt_thcd

ROOT = Path(__file__).parent
LOG_PATH = ROOT / "LOG.md"
TC_SUPPORTED = {
    "Al", "Co", "Cr", "Fe", "Ni", "Mn", "Mo", "Nb", "Ta", "Ti",
    "V", "W", "Zr", "Hf", "Cu", "Si", "B", "C", "N", "Y", "Re", "Ru",
}


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module {name} from {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_props_mod = _load_module(ROOT / "1_find_stoich_props_parallel_for_Class.py", "cafer_props")
compute_basic_features = _props_mod.compute_basic_features
predict_ys_1300c = _props_mod.predict_ys_1300c


def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with LOG_PATH.open("a") as f:
        f.write(f"- [{ts}] {msg}\n")


def parse_formula(formula: str) -> dict:
    """Parse formula like Al0.3NbTaTi1.4Zr1.3 into element->amount."""
    if not isinstance(formula, str):
        return {}
    tokens = re.findall(r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)", formula)
    comp = {}
    for el, amt in tokens:
        val = float(amt) if amt not in ("", None) else 1.0
        comp[el] = comp.get(el, 0.0) + val
    return comp


def normalize_row(row: pd.Series, element_cols: list[str]) -> pd.Series:
    comp = row[element_cols].astype(float).fillna(0.0)
    total = comp.sum()
    if total > 0:
        comp = comp / total
    row[element_cols] = comp
    return row


def build_from_formula(df: pd.DataFrame, col: str) -> tuple[pd.DataFrame, list[str]]:
    parsed = [parse_formula(v) for v in df[col].fillna("")]
    elements = sorted({el for comp in parsed for el in comp})
    comp_df = pd.DataFrame([{el: comp.get(el, 0.0) for el in elements} for comp in parsed])
    comp_df = comp_df.apply(lambda r: normalize_row(r, elements), axis=1)
    return comp_df, elements


def featurize(df: pd.DataFrame, element_cols: list[str], tag: str, do_tc: bool) -> pd.DataFrame:
    log(f"Featurizing {tag}: {len(df)} rows, {len(element_cols)} element columns")
    work = df.copy()
    work[element_cols] = work[element_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    work = work.apply(lambda r: normalize_row(r, element_cols), axis=1)

    work = compute_basic_features(work, element_cols, log_fn=log)
    work = predict_ys_1300c(work, element_cols, log_fn=log)

    def is_tc_supported(row):
        active = [el for el in element_cols if row[el] > 0]
        return set(active).issubset(TC_SUPPORTED)

    if do_tc:
        mask = work.apply(is_tc_supported, axis=1)
        tc_rows = work[mask].copy()
        no_tc = work[~mask].copy()
        if len(tc_rows) > 0:
            tc_rows = compute_solidus_liquidus(tc_rows, element_cols, log_fn=log)
            tc_rows = compute_rt_thcd(tc_rows, element_cols, log_fn=log)
        else:
            log(f"No TC-supported rows for {tag}")
    else:
        tc_rows = work.iloc[0:0].copy()
        no_tc = work.copy()

    # Ensure TC columns exist
    for col in ["Solidus_K_TC", "Liquidus_K_TC", "Solidus_C_TC", "Liquidus_C_TC", "THCD_25C_W_mK"]:
        if col not in tc_rows.columns:
            tc_rows[col] = np.nan
        if col not in no_tc.columns:
            no_tc[col] = np.nan

    work = pd.concat([tc_rows, no_tc], axis=0).sort_index()

    return work


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", default="borg2019,codhem,lietal", help="Comma list: borg2019,codhem,lietal")
    parser.add_argument("--skip-tc", action="store_true", help="Skip TC solidus/liquidus and RT THCD")
    args = parser.parse_args()

    do_tc = not args.skip_tc
    selected = {s.strip() for s in args.datasets.split(",") if s.strip()}

    log(f"Starting featurization run (datasets={sorted(selected)}, do_tc={do_tc})")

    # borg2019.csv
    if "borg2019" in selected:
        borg_path = ROOT / "borg2019.csv"
        borg = pd.read_csv(borg_path)
        borg_comp, borg_elements = build_from_formula(borg, "FORMULA")
        borg_feat = featurize(pd.concat([borg, borg_comp], axis=1), borg_elements, "borg2019.csv", do_tc)
        borg_out = ROOT / "borg2019_featurized.csv"
        borg_feat.to_csv(borg_out, index=False)
        log(f"Wrote {borg_out}")

    # codhem.csv
    if "codhem" in selected:
        cod_path = ROOT / "codhem.csv"
        cod = pd.read_csv(cod_path)
        cod_comp, cod_elements = build_from_formula(cod, "Composition")
        cod_feat = featurize(pd.concat([cod, cod_comp], axis=1), cod_elements, "codhem.csv", do_tc)
        cod_out = ROOT / "codhem_featurized.csv"
        cod_feat.to_csv(cod_out, index=False)
        log(f"Wrote {cod_out}")

    # lietal.xlsx (element columns already present)
    if "lietal" in selected:
        li_path = ROOT / "lietal.xlsx"
        li = pd.read_excel(li_path)
        start_idx = li.columns.get_loc("Ag")
        li_elements = list(li.columns[start_idx:])
        li_feat = featurize(li, li_elements, "lietal.xlsx", do_tc)
        li_out = ROOT / "lietal_featurized.xlsx"
        li_feat.to_excel(li_out, index=False)
        log(f"Wrote {li_out}")

    log("Featurization run complete")


if __name__ == "__main__":
    main()
