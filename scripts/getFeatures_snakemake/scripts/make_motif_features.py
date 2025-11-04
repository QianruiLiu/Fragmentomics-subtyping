#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Summarize end-motif (k-mer) composition from either a wide table
#        (many k-mer columns per row) or a long table (one k-mer column + value).
#        Outputs a single-row TSV with entropy and grouped prefix fractions.
# Usage:
#   python make_motif_features.py input.tsv output.tsv \
#     --kmer_regex "(?i)^kmer$" --value_col count \
#     --prefixes_CG CG GC --prefixes_AT AA TT AT TA --prefixes_CCGG CC GG
# Notes:
#   * Auto-detects wide vs long format by scanning for k-mer-like column names.
#   * In wide mode, optionally weights rows by a "count" column if present.
#   * In long mode, aggregates values by k-mer and normalizes to probabilities.

import sys, argparse, re, numpy as np, pandas as pd

def looks_like_kmer(s):
    # Accept 2–6 bp DNA tokens (A/C/G/T/N), case-insensitive
    return bool(re.fullmatch(r"[ACGTNacgtn]{2,6}", str(s)))

ap = argparse.ArgumentParser()
ap.add_argument("tsv_in")
ap.add_argument("out_tsv")
ap.add_argument("--kmer_regex", type=str, default=r"(?i)^kmer$")
ap.add_argument("--value_col", type=str, default="count")
ap.add_argument("--prefixes_CG", nargs="*", default=["CG","GC"])
ap.add_argument("--prefixes_AT", nargs="*", default=["AA","TT","AT","TA"])
ap.add_argument("--prefixes_CCGG", nargs="*", default=["CC","GG"])
args = ap.parse_args()

t = pd.read_csv(args.tsv_in, sep="\t", engine="python")

# ---- Two paths: wide vs long ----
kmer_cols_wide = [c for c in t.columns if looks_like_kmer(c)]
# Heuristic: count-like number of k-mer columns (e.g., ~256 for 4-mers); keep threshold lenient
is_wide = len(kmer_cols_wide) >= 16

if is_wide:
    # Wide format: each row is an interval, many k-mer columns hold frequencies or counts
    kmers = [str(c).upper() for c in kmer_cols_wide]
    T = t.copy()
    T_cols = dict(zip(kmer_cols_wide, kmers))
    T.rename(columns=T_cols, inplace=True)

    if "count" in T.columns:
        # Weighted sum across rows: sum_row(freq_or_count * row_count)
        weights = T["count"].astype(float).fillna(0.0).to_numpy()
        mat = T[kmers].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()
        # Broadcast per-row weights
        weighted = (mat.T * weights).T
        v = weighted.sum(axis=0)
    else:
        # No explicit count column: simple sum across rows
        v = T[kmers].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=0).to_numpy()

    kmer_list = kmers
    total = float(np.nansum(v))
    p = v / total if total > 0 else v.astype(float)

else:
    # Long format: must have a k-mer column and a value (count/freq) column
    kmer_col = None
    for c in t.columns:
        if re.match(args.kmer_regex, str(c)):
            kmer_col = c; break
    if kmer_col is None:
        for c in t.columns:
            if 'kmer' in str(c).lower():
                kmer_col = c; break
    if kmer_col is None:
        sys.exit("ERROR: k-mer column not found; wide table with k-mer columns was not detected. Check input or adjust --kmer_regex.")

    val_col = args.value_col if args.value_col in t.columns else None
    if val_col is None:
        for cand in ["count","counts","freq","frequency","value","n","total","sum"]:
            if cand in t.columns: val_col=cand; break
    if val_col is None:
        sys.exit("ERROR: value column (count/freq) not found. Specify it via --value_col.")

    # Aggregate by k-mer and normalize later
    g = t.groupby(kmer_col)[val_col].sum(numeric_only=True).reset_index()
    g.columns = ["kmer","value"]
    g["kmer"] = g["kmer"].astype(str).str.upper()
    v = g["value"].astype(float).to_numpy()
    kmer_list = g["kmer"].tolist()
    total = float(np.nansum(v))
    p = v / total if total > 0 else v

# ---- Features: entropy + grouped prefix fractions ----
def frac(prefixes):
    # Fraction of total probability mass for k-mers starting with any of the prefixes
    mask = np.array([any(k.startswith(px.upper()) for px in prefixes) for k in kmer_list])
    return float(np.nansum(p[mask])) if total>0 else np.nan

# Shannon entropy (base-2), guard against log(0)
entropy = float(-np.nansum(p*np.log2(np.clip(p,1e-12,1)))) if total>0 else np.nan

pd.DataFrame([{
  "motif_entropy": entropy,
  "frac_CGorGC": frac(args.prefixes_CG),
  "frac_AT_rich": frac(args.prefixes_AT),
  "frac_CC_GG": frac(args.prefixes_CCGG)
}]).to_csv(args.out_tsv, sep="\t", index=False)
