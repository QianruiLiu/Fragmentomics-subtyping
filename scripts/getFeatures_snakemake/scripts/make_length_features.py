#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Compute fragment-length features from either a histogram-style table
#        (multiple columns like "100_150", "160_220", ...) or a summary table
#        (columns like mean/median/stdev/count and optional s150). Writes a
#        single-row TSV with derived features.
# Usage:
#   python make_length_features.py input.bed output.tsv \
#     --short 100 150 --mono 160 220 --di 300 400 --short_threshold 150
# Notes:
#   * Auto-detects input format by scanning for length-bin columns.
#   * In histogram mode, features are fractions over user-specified bins.
#   * In summary mode, features are weighted aggregates over intervals.

import sys, argparse, re, numpy as np, pandas as pd

ap = argparse.ArgumentParser(description="Compute fragment-length features (histogram or summary format).")
ap.add_argument("bed_in")
ap.add_argument("out_tsv")
# Binning for histogram mode (inclusive)
ap.add_argument("--short", nargs=2, type=int, default=[100,150])
ap.add_argument("--mono",  nargs=2, type=int, default=[160,220])
ap.add_argument("--di",    nargs=2, type=int, default=[300,400])
# Summary-mode hint (informational only; does not change how s150 is read)
ap.add_argument("--short_threshold", type=int, default=150, help="threshold for s150 meaning (<=bp)")
args = ap.parse_args()

df = pd.read_csv(args.bed_in, sep="\t")

# Normalize column names: strip whitespace and lower-case for robust matching
df.columns = [str(c).strip() for c in df.columns]
lower_map = {c: c.lower() for c in df.columns}
df.rename(columns=lower_map, inplace=True)

def is_len_bin(name:str)->bool:
    # Detect columns that look like length bins: "100_150", "100-150", or "len_100_150"
    s=str(name)
    return bool(re.search(r'(\d+)[_-](\d+)', s)) or bool(re.search(r'len[_-]?(\d+)[_-](\d+)', s, re.I))

# ------ Decide branch: histogram vs summary format ------
bin_cols = [c for c in df.columns if is_len_bin(c)]

if len(bin_cols) >= 3:
    # ========= Histogram mode =========
    # Convert selected bin columns to numeric matrix (rows = intervals/entries, cols = bins)
    arr = df[bin_cols].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)

    # Parse each bin's lower/upper bounds from the column name
    bounds=[]
    for c in bin_cols:
        m = re.search(r'(\d+)[_-](\d+)', str(c))
        if m: bounds.append((int(m.group(1)), int(m.group(2))))
        else: bounds.append((None,None))

    def pick(lo,hi):
        # Sum counts across all bins fully inside [lo, hi]
        idx=[i for i,(a,b) in enumerate(bounds) if a is not None and a>=lo and b<=hi]
        if len(idx)==0:
            return np.zeros(arr.shape[0], dtype=float)
        return np.nansum(arr[:, idx], axis=1)

    S = pick(*args.short); M = pick(*args.mono); D = pick(*args.di)
    T = np.nansum(arr, axis=1)

    S_sum, M_sum, D_sum, T_sum = map(lambda x: float(np.nansum(x)), [S,M,D,T])

    short_frac = S_sum / T_sum if T_sum>0 else np.nan
    mono_frac = M_sum / T_sum if T_sum>0 else np.nan
    di_frac = D_sum / T_sum if T_sum>0 else np.nan
    short_over_mono = S_sum / M_sum if M_sum>0 else np.nan

    # Global length-distribution entropy over all bins
    p = np.nansum(arr, axis=0)
    p = p / np.nansum(p) if np.nansum(p)>0 else p
    length_entropy = float(-np.nansum(p*np.log2(np.clip(p,1e-12,1)))) if np.nansum(p)>0 else np.nan

    out = pd.DataFrame([{
        "short_frac": short_frac,
        "mono_frac": mono_frac,
        "di_frac": di_frac,
        "short_over_mono": short_over_mono,
        "length_entropy": length_entropy,
        "mode": "histogram"
    }])

else:
    # ========= Summary mode =========
    # Support common aliases for required summary columns
    col_aliases = {
        "mean": ["mean", "avg", "average"],
        "median": ["median"],
        "stdev": ["stdev", "stddev", "std"],
        "count": ["count", "n", "total"]
    }
    resolved = {}
    for key, cands in col_aliases.items():
        for c in cands:
            if c in df.columns:
                resolved[key] = c
                break
        if key not in resolved:
            sys.exit(f"ERROR: summary mode requires columns {list(col_aliases.keys())} (missing {key})")

    # Optional s150 column: accept exact 's150' or relaxed variants (e.g., 's_150')
    s150_col = None
    if "s150" in df.columns:
        s150_col = "s150"
    else:
        for c in df.columns:
            if re.fullmatch(r"s[\W_]*150", c):
                s150_col = c
                break

    # Coerce core summary columns to numeric
    for c in [resolved["mean"], resolved["median"], resolved["stdev"], resolved["count"]]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    d = df.copy()
    d = d[(d[resolved["count"]]>0)]
    # Treat negative min/max as missing if present
    if "min" in d.columns:
        d["min"] = pd.to_numeric(d["min"], errors="coerce")
    if "max" in d.columns:
        d["max"] = pd.to_numeric(d["max"], errors="coerce")

    total = float(d[resolved["count"]].sum())
    if total<=0:
        out = pd.DataFrame([{
            "weighted_mean_len": np.nan,
            "weighted_median_len": np.nan,
            "pooled_sd_len": np.nan,
            "short_frac_s150": (np.nan if s150_col else None),
            "total_fragments": 0,
            "n_intervals": int(d.shape[0]),
            "global_min_len": (float("nan") if "min" in d.columns else None),
            "global_max_len": (float("nan") if "max" in d.columns else None),
            "mode": "summary"
        }])
    else:
        # Weighted mean by counts
        w_mean = float((d[resolved["mean"]]*d[resolved["count"]]).sum()/total)

        # Approximate weighted median: repeat medians by their counts (can be large)
        med = d[[resolved["median"], resolved["count"]]].dropna()
        weighted_median = np.nan
        try:
            rep = np.repeat(med[resolved["median"]].to_numpy(),
                            med[resolved["count"]].astype(int).to_numpy())
            weighted_median = float(np.median(rep)) if rep.size>0 else np.nan
        except Exception:
            pass

        # Pooled standard deviation (within + between)
        var_within = ((d[resolved["stdev"]]**2) * (d[resolved["count"]]-1)).sum()
        var_between = (d[resolved["count"]] * (d[resolved["mean"]] - w_mean)**2).sum()
        pooled_var = (var_within + var_between) / max(total-1, 1)
        pooled_sd = float(np.sqrt(pooled_var))

        # Optional s150 weighted fraction if column exists
        short_frac = np.nan
        if s150_col:
            d[s150_col] = pd.to_numeric(d[s150_col], errors="coerce")
            short_frac = float((d[s150_col]*d[resolved["count"]]).sum()/total)

        gmin = float(d["min"][d["min"]>=0].min()) if "min" in d.columns else np.nan
        gmax = float(d["max"][d["max"]>=0].max()) if "max" in d.columns else np.nan

        out = pd.DataFrame([{
            "weighted_mean_len": w_mean,
            "weighted_median_len": weighted_median,
            "pooled_sd_len": pooled_sd,
            "short_frac_s150": short_frac,
            "total_fragments": int(total),
            "n_intervals": int(d.shape[0]),
            "global_min_len": gmin,
            "global_max_len": gmax,
            "mode": "summary"
        }])

# If you need a fixed schema, you could add missing histogram columns as NA here (optional).
out.to_csv(args.out_tsv, sep="\t", index=False)
