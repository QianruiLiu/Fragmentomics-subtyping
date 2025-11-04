#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Merge per-sample feature tables produced under results/<sample>/<signal>/<site_set>/
#        into one wide matrix (rows = samples, columns = features). Supports multiple signals
#        (e.g., cleavage_profile, fraglen, motif, wps) and two site-sets (AD/NE). Optionally
#        joins a labels TSV on sample_id. When duplicate rows for a sample exist, prefer
#        non-zero values per column.
# Usage:
#   python merge_features.py --results /path/to/results \
#     --signals cleavage_profile fraglen motif wps \
#     --out merged/features.tsv [--labels labels.tsv]
# Notes:
#   * Expects each leaf features file to be named *.features.tsv.
#   * Column names in each features file will be prefixed as "{site_set}_{signal}_".
#   * Site-sets are assumed to be "AD" and "NE" subfolders under each signal.

import argparse, os, re, glob
import pandas as pd

def base_sample_id(name: str) -> str:
    # Normalize by trimming a trailing status tag from the directory name: supports old (ARPC/NEPC) and new (AD/NE)
    return re.sub(r'_(ARPC|NEPC|AD|NE)$', '', name, flags=re.IGNORECASE)

def read_features_dir(sample_dir: str, signal: str, site_set: str):
    """
    Read *.features.tsv for a given sample, signal, and site_set (one of {'AD','NE'}).
    Returns a single-row DataFrame with columns prefixed as "{site_set}_{signal}_".
    """
    subdir = os.path.join(sample_dir, signal, site_set)
    if not os.path.isdir(subdir):
        return None

    # Allow one or multiple features files; if multiple, merge side-by-side (usually there's only one)
    tsvs = sorted(glob.glob(os.path.join(subdir, "*.features.tsv")))
    if not tsvs:
        return None

    dfs = []
    for tsv in tsvs:
        try:
            df = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue

        # Common case: features.tsv is "one row, many columns"; if a sample_id column exists, keep the first row
        if df.shape[0] > 1:
            # If 'sample_id' is present, group and take the first row; otherwise just take the first row
            if "sample_id" in df.columns:
                df = df.groupby("sample_id").head(1).reset_index(drop=True)
            else:
                df = df.head(1)

        # Drop id-like columns (we add sample_id externally at the row level)
        drop_cols = [c for c in df.columns if c.lower() in {"sample","sample_id","id","label"}]
        df = df.drop(columns=drop_cols, errors="ignore")

        # Add a prefix such as AD_cleavage_profile_center_over_shoulder
        df = df.add_prefix(f"{site_set}_{signal}_")
        dfs.append(df)

    if not dfs:
        return None

    # Merge all features in this directory horizontally (most often just one file)
    out = pd.concat(dfs, axis=1)
    return out

def collect_one_sample(sample_dir: str, signals):
    """
    Build the aggregated single row for a sample, spanning both site-sets (AD/NE) across all signals.
    """
    all_parts = []
    for site_set in ("AD","NE"):
        for signal in signals:
            df = read_features_dir(sample_dir, signal, site_set)
            if df is not None and not df.empty:
                all_parts.append(df)

    if not all_parts:
        return None

    row = pd.concat(all_parts, axis=1)
    return row

# -------- Prefer-nonzero reducer used when collapsing duplicate sample rows --------
def _prefer_nonzero(series: pd.Series):
    """
    Column-wise reducer when merging multiple rows per sample_id:
    - Prefer a non-zero (and non-NA, non-empty) value.
    - If all are zero or non-numeric, fall back to the first non-empty element.
    """
    s = series.dropna()
    s = s[s.astype(str) != ""]
    if s.empty:
        return pd.NA
    # Try numeric evaluation to identify non-zero entries
    sn = pd.to_numeric(s, errors="coerce")
    nz = sn[(~sn.isna()) & (sn != 0)]
    if len(nz):
        idx = nz.index[0]
        return s.loc[idx]
    # No non-zero values; return the first non-empty
    return s.iloc[0]
# ----------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Merge AD/NE features (multi-signal) into a sample×feature matrix.")
    ap.add_argument("--results", required=True, help="Results root directory, e.g., /home/student/getFeatures_snakemake/results")
    ap.add_argument("--signals", nargs="+", default=["cleavage_profile","fraglen","motif","wps"],
                    help="Names of signal subdirectories to merge (edit as needed).")
    ap.add_argument("--labels", default=None, help="Optional labels TSV containing 'sample_id' plus covariates/targets.")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    args = ap.parse_args()

    # Discover per-sample folders under results/
    sample_dirs = [os.path.join(args.results, d) for d in os.listdir(args.results)
                   if os.path.isdir(os.path.join(args.results, d))]
    records = []
    index_names = []

    for sdir in sorted(sample_dirs):
        sname = os.path.basename(sdir)
        sid = base_sample_id(sname)
        row = collect_one_sample(sdir, args.signals)
        if row is None or row.empty:
            continue
        row.insert(0, "sample_id", sid)
        records.append(row)
        index_names.append(sid)

    if not records:
        raise SystemExit("No features collected. Check --results and folder layout.")

    features_df = pd.concat(records, axis=0, ignore_index=True)

    # Collapse potential duplicate rows per sample_id using the prefer-nonzero reducer
    features_df = (features_df.groupby("sample_id", as_index=False)
                              .agg(_prefer_nonzero))

    # Optional: join labels if provided
    if args.labels:
        lab = pd.read_csv(args.labels, sep="\t")
        # Keep 'sample_id' as the key; other columns (label, tumor_fraction, batch, etc.) are carried along
        if "sample_id" not in lab.columns:
            raise SystemExit("Labels file must contain 'sample_id' column.")
        features_df = features_df.merge(lab, on="sample_id", how="left")

    # Optional tidy column order: sample_id first, then AD_* columns, then NE_*, then any others
    cols = list(features_df.columns)
    cols.remove("sample_id")
    cols_sorted = (["sample_id"] +
                   sorted([c for c in cols if c.startswith("AD_")]) +
                   sorted([c for c in cols if c.startswith("NE_")]) +
                   [c for c in cols if not (c.startswith("AD_") or c.startswith("NE_"))])
    features_df = features_df.reindex(columns=cols_sorted)

    # Write output
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    features_df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote {args.out} with shape {features_df.shape}")

if __name__ == "__main__":
    main()
