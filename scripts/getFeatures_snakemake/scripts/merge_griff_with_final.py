#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Normalize sample names (with flexible LuCaP patterns) and merge two TSVs
#        (finaletoolkit = left/base, griffin = right) by the normalized key.
#        Keeps columns from A, adds only non-duplicate columns from B, and writes
#        a merge report listing unmatched keys on each side.
# Usage:
#   python merge_two_feature_tables.py \
#     --finaletoolkit A.tsv --griffin B.tsv --out merged.tsv
# Notes:
#   * Sample column is auto-detected (sample/samples/sample_id/sampleid/id; else first column).
#   * Name normalization supports both "LuCaP_<num>..." and "<num>..._LuCaP" with optional CR, _sub, _repN.
#   * Duplicate rows per key are de-duplicated by keeping the last occurrence in each table.

import re
import argparse
import pandas as pd
from pathlib import Path

# --------------------------
# 1) Identify the sample column
# --------------------------
def detect_sample_col(df: pd.DataFrame) -> str:
    lower_map = {str(c).lower(): c for c in df.columns}
    for key in ["sample", "samples", "sample_id", "sampleid", "id"]:
        if key in lower_map:
            return lower_map[key]
    # Otherwise, use the first column
    return df.columns[0]

# --------------------------
# 2) Normalize sample names
# Unify both orders “LuCaP_###…” and “###…_LuCaP”.
# Support CR, sub-index (e.g., _1/_2), and replication tags (repN).
# --------------------------
_COMMON_SUFFIXES = [
    ".tsv", ".csv", ".txt",
    ".bam", ".cram", ".bw", ".bed",
    ".gz",
    "_recal", ".recal",
    "_filtered", ".filtered",
    "_dedup", ".dedup"
]

def strip_common_suffixes(s: str) -> str:
    out = s
    changed = True
    while changed:
        changed = False
        for suf in _COMMON_SUFFIXES:
            if out.endswith(suf):
                out = out[: -len(suf)]
                changed = True
    return out

def _cleanup_base(s: str) -> str:
    """Remove path/extension, normalize separators, drop common technical suffixes."""
    s = str(s)
    s = Path(s).name  # drop directory path
    s = s.strip().strip('\'"')
    s = strip_common_suffixes(s)
    s = s.replace("-", "_")
    s = re.sub(r"[ \t]+", "_", s)
    s = re.sub(r"__+", "_", s)
    s = s.strip("_")
    return s

def normalize_name(x: str) -> str:
    if pd.isna(x):
        return x
    s = _cleanup_base(x)

    # Use lowercase only for matching; keep canonical capitalization ("LuCaP") in output
    s_low = s.lower()

    # Allowed optional fragments:
    #  - num: digits such as 136, 35, 176
    #  - cr : optional "CR" immediately after digits (e.g., 35CR, 176CR)
    #  - sub: optional sub-index like _1, _2, _4
    #  - rep: optional replicate like _rep2
    # Two possible orders:
    #  A) LuCaP_<num>(CR)?(_sub)?(_rep\d+)?
    #  B) <num>(CR)?(_sub)?(_rep\d+)?_LuCaP
    # Note: In real data, sub usually precedes rep (e.g., 173_1_rep2). We match that first.
    fwd = re.compile(r'^lucap_(?P<num>\d+)(?P<cr>cr)?(?P<sub>_\d+)?(?P<rep>_rep\d+)?$', re.IGNORECASE)
    rev = re.compile(r'^(?P<num>\d+)(?P<cr>cr)?(?P<sub>_\d+)?(?P<rep>_rep\d+)?_lucap$', re.IGNORECASE)

    m = fwd.match(s_low)
    if not m:
        m = rev.match(s_low)

    # Compatibility: if the rare order rep-before-sub appears, try an alternate pattern
    if not m:
        fwd_alt = re.compile(r'^lucap_(?P<num>\d+)(?P<cr>cr)?(?P<rep>_rep\d+)?(?P<sub>_\d+)?$', re.IGNORECASE)
        rev_alt = re.compile(r'^(?P<num>\d+)(?P<cr>cr)?(?P<rep>_rep\d+)?(?P<sub>_\d+)?_lucap$', re.IGNORECASE)
        m = fwd_alt.match(s_low) or rev_alt.match(s_low)

    if m:
        num = m.group("num")
        cr = m.group("cr") or ""
        sub = m.group("sub") or ""
        rep = m.group("rep") or ""
        # Canonical form: LuCaP_<num>[CR][_sub][_repN]
        cr_up = "CR" if cr else ""
        return f"LuCaP_{num}{cr_up}{sub}{rep}"

    # Fallback: already starts with LuCaP but didn't match strict pattern (e.g., LuCaP_170_2). Return lightly cleaned.
    if s_low.startswith("lucap_"):
        return "LuCaP_" + s[6:]

    # Fallback: ends with _LuCaP but didn't match strict pattern
    if s_low.endswith("_lucap"):
        core = s[:-6].strip("_")
        return f"LuCaP_{core}"

    # Last resort: return the cleaned string unchanged
    return s

# --------------------------
# 3) Main routine
# --------------------------
def main(final_path: str, griffin_path: str, out_path: str):
    dfA = pd.read_csv(final_path, sep="\t", dtype=str)
    dfB = pd.read_csv(griffin_path, sep="\t", dtype=str)

    # Detect sample column in each table
    a_sample = detect_sample_col(dfA)
    b_sample = detect_sample_col(dfB)

    # Normalize sample names into a join key "_key"
    dfA["_key"] = dfA[a_sample].map(normalize_name)
    dfB["_key"] = dfB[b_sample].map(normalize_name)

    # Drop duplicates per key, keeping the last occurrence
    dfA = dfA.dropna(subset=["_key"]).drop_duplicates(subset=["_key"], keep="last")
    dfB = dfB.dropna(subset=["_key"]).drop_duplicates(subset=["_key"], keep="last")

    # Choose columns to add from B: skip columns already present in A (except the sample columns and _key)
    a_cols = set(dfA.columns)
    b_cols_to_add = []
    for col in dfB.columns:
        if col in {"_key", b_sample}:
            continue
        if col in a_cols:
            # Per requirement: do not add suffixes or duplicate-named columns; skip if name already exists in A
            continue
        b_cols_to_add.append(col)

    dfB_add = dfB[["_key"] + b_cols_to_add].copy()

    # Left join with A as the base table
    merged = pd.merge(dfA, dfB_add, on="_key", how="left", copy=False)

    # Clean up helper column: keep A's sample column(s); drop the join key
    merged = merged.drop(columns=["_key"])

    # Write merged matrix
    merged.to_csv(out_path, sep="\t", index=False)

    # Produce a match report (lists keys unique to each side)
    base = Path(out_path).with_suffix("")
    report_path = f"{base}_merge_report.txt"
    a_keys = set(dfA["_key"])
    b_keys = set(dfB["_key"])
    only_in_a = sorted(a_keys - b_keys)
    only_in_b = sorted(b_keys - a_keys)

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(f"Total in A (finaletoolkit): {len(a_keys)}\n")
        f.write(f"Total in B (griffin)     : {len(b_keys)}\n")
        f.write(f"Matched                  : {len(a_keys & b_keys)}\n\n")
        f.write("Only in A (no match in B):\n")
        for k in only_in_a:
            f.write(f"  {k}\n")
        f.write("\nOnly in B (not in A base):\n")
        for k in only_in_b:
            f.write(f"  {k}\n")

    print(f"[OK] Merged → {out_path}")
    print(f"[OK] Report → {report_path}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Merge finaletoolkit & griffin feature tables by normalized sample names (handles 'LuCaP' prefix/suffix).")
    ap.add_argument("--finaletoolkit", required=True, help="TSV file from finaletoolkit (base table).")
    ap.add_argument("--griffin", required=True, help="TSV file from griffin.")
    ap.add_argument("--out", required=True, help="Output TSV for merged feature matrix.")
    args = ap.parse_args()
    main(args.finaletoolkit, args.griffin, args.out)
