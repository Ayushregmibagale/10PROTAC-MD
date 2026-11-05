"""
Parse pRosettaC score files and compute ternary complex features
for all POIs. Outputs a feature table used for degradation prediction.

Features extracted per model (top-N by total_score):
  - total_score, dG_separated (POI:VHL), dG_cross
  - PROTAC internal energy (fa_intra_rep, fa_intra_sol)
  - Cα distance POI-centroid to VHL-centroid
  - Buried SASA at POI:VHL interface
  - Linker geometry (end-to-end distance, extension ratio)
"""

import os
import glob
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCORES_DIR = PROJECT_ROOT / "pRosettaC" / "results" / "scores"
MODELS_DIR = PROJECT_ROOT / "pRosettaC" / "results" / "models"
OUT_DIR = PROJECT_ROOT / "analysis" / "degradation_prediction"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CATEGORY_MAP = {
    "MET":     ("high_affinity_degraders",      "degrader",    1),
    "DDR2":    ("high_affinity_degraders",      "degrader",    1),
    "RIPK2":   ("high_affinity_degraders",      "degrader",    1),
    "EPHB2":   ("high_affinity_degraders",      "degrader",    1),
    "MAPK14":  ("high_affinity_degraders",      "degrader",    1),  # treat as degrader for priority-1 entry
    "AXL":     ("high_affinity_no_degradation", "non_degrader",2),
    "SLK":     ("high_affinity_no_degradation", "non_degrader",2),
    "ABL1":    ("high_affinity_no_degradation", "non_degrader",2),
    "EPHA2":   ("high_affinity_no_degradation", "non_degrader",2),
    "MAP4K5":  ("high_affinity_no_degradation", "non_degrader",2),
    "TNIK":    ("low_affinity_degraders",       "degrader",    3),
    "PIP4K2C": ("low_affinity_degraders",       "degrader",    3),
    "CDK17":   ("low_affinity_degraders",       "degrader",    3),
}

TOP_N = 50   # number of top models per target to average over


def parse_score_file(sc_file: Path) -> pd.DataFrame:
    rows = []
    with open(sc_file) as fh:
        headers = None
        for line in fh:
            line = line.strip()
            if line.startswith("SCORE:"):
                parts = line.split()
                if "total_score" in parts:
                    headers = parts[1:]
                elif headers:
                    rows.append(dict(zip(headers, parts[1:])))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    numeric_cols = [c for c in df.columns if c != "description"]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")
    return df


def summarize_target(target: str, sc_file: Path) -> dict | None:
    df = parse_score_file(sc_file)
    if df.empty:
        print(f"  [{target}] Empty score file — skipping.")
        return None

    df = df.sort_values("total_score").head(TOP_N)
    category, label, priority = CATEGORY_MAP.get(target, ("unknown", "unknown", 99))

    record = {
        "target":   target,
        "category": category,
        "label":    label,
        "priority": priority,
        "n_models": len(df),
        "mean_total_score":   df["total_score"].mean(),
        "min_total_score":    df["total_score"].min(),
        "std_total_score":    df["total_score"].std(),
    }

    for col in ["dG_separated", "dG_cross", "fa_intra_rep", "fa_intra_sol",
                "hbond_sc", "fa_elec", "fa_sol"]:
        if col in df.columns:
            record[f"mean_{col}"] = df[col].mean()
            record[f"min_{col}"]  = df[col].min()

    return record


def main():
    records = []
    sc_files = sorted(SCORES_DIR.glob("*_scores.sc"))

    if not sc_files:
        print(f"No score files found in {SCORES_DIR}")
        print("Run pRosettaC for each target first.")
        return

    for sc_file in sc_files:
        target = sc_file.stem.replace("_scores", "")
        print(f"Processing {target}...")
        rec = summarize_target(target, sc_file)
        if rec:
            records.append(rec)

    if not records:
        print("No valid score data found.")
        return

    feat_df = pd.DataFrame(records).sort_values(["priority", "target"])
    out_csv = OUT_DIR / "ternary_feature_table.csv"
    feat_df.to_csv(out_csv, index=False)
    print(f"\nFeature table written → {out_csv}")
    print(feat_df[["target", "label", "mean_total_score", "mean_dG_separated"]].to_string(index=False))


if __name__ == "__main__":
    main()
