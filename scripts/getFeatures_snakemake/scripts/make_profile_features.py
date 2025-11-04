#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Convert a tidy meta-profile (deepTools-like output) into (1) a smoothed profile
#        and (2) engineered features: central/shoulder means, contrast, center-over-shoulder,
#        amplitude in a core window, and phasing strength from four flank segments.
# Usage:
#   python make_profile_features.py input.tidy.tsv out_prefix \
#     --binsize 15 --center -30 30 --shoulders -150 -75 75 150 \
#     --amp_core -150 150 --phase -220 -170 -170 -120 120 170 170 220 --smooth 165
# Notes:
#   * Accepts two common deepTools exports: (A) long/tidy with 'x' column, (B) matrix-like row 2.. as y.
#   * Smoothing uses Savitzky–Golay with a bp-sized window converted to bins (0 = off).
#   * Positions are in bp; masks are half-open [lo, hi).

import sys, math, argparse
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

ap = argparse.ArgumentParser()
ap.add_argument("tidy_tsv")
ap.add_argument("out_prefix")
ap.add_argument("--binsize", type=int, required=True, help="Bin size in bp for the meta-profile.")
ap.add_argument("--center", nargs=2, type=int, default=[-30,30], help="Center window [lo hi) for central mean.")
ap.add_argument("--shoulders", nargs=4, type=int, default=[-150,-75,75,150], help="Shoulder bands: L1 L2 R1 R2 (two half-open ranges).")
ap.add_argument("--amp_core", nargs=2, type=int, default=[-150,150], help="Core window for amplitude (max-min).")
ap.add_argument("--phase", nargs=8, type=int, default=[-220,-170,-170,-120,120,170,170,220],
                help="Phasing segments: Lp1 Lp2 Lv1 Lv2 Rv1 Rv2 Rp1 Rp2 (four half-open ranges).")
ap.add_argument("--smooth", type=int, default=0, help="Savitzky–Golay window (bp), 0=off.")
args = ap.parse_args()

df = pd.read_csv(args.tidy_tsv, sep="\t")

# Support two deepTools export styles:
# (A) Has 'x' column -> use first non-x column as the averaged profile
# (B) Matrix-like table -> reconstruct x from binsize and array length
if "x" in df.columns and df.shape[1] >= 2:
    x = df["x"].to_numpy()
    y = df[df.columns[df.columns != "x"][0]].to_numpy(dtype=float)
    pos = x
else:
    y = df.iloc[1, 2:].astype(float).to_numpy()
    n = y.size
    half = (n * args.binsize) // 2
    pos = np.arange(-half, -half + n * args.binsize, args.binsize)

# Optional smoothing (Savitzky–Golay). Require window >= 3 bins and odd.
if args.smooth and args.smooth > args.binsize:
    win_bins = max(3, int(round(args.smooth / args.binsize)) | 1)  # force odd
    y = savgol_filter(y, window_length=win_bins, polyorder=3, mode="interp")

def mask(lo, hi):
    # Half-open interval [lo, hi)
    return (pos >= lo) & (pos < hi)

def mean_safe(v):
    # Mean that ignores NaNs and returns NaN if empty
    v = v[np.isfinite(v)]
    return float("nan") if v.size == 0 else float(np.mean(v))

# Windows and feature blocks
center = mask(*args.center)
L1, L2, R1, R2 = args.shoulders
shoulder = mask(L1, L2) | mask(R1, R2)
amp_core = mask(*args.amp_core)

central_mean = mean_safe(y[center])
shoulder_mean = mean_safe(y[shoulder])
contrast = central_mean - shoulder_mean if np.isfinite([central_mean, shoulder_mean]).all() else float("nan")
ratio = central_mean / shoulder_mean if (np.isfinite([central_mean, shoulder_mean]).all() and not math.isclose(shoulder_mean, 0.0)) else float("nan")

# Amplitude within the core window (max - min)
core_vals = y[amp_core]; core_vals = core_vals[np.isfinite(core_vals)]
amplitude = float("nan") if core_vals.size == 0 else float(np.max(core_vals) - np.min(core_vals))

# Phasing strength: average (peak - valley) on left and right
Lp1, Lp2, Lv1, Lv2, Rv1, Rv2, Rp1, Rp2 = args.phase
def seg(lo, hi):
    vv = y[mask(lo, hi)]; vv = vv[np.isfinite(vv)]
    return float("nan") if vv.size == 0 else float(np.mean(vv))
left_peak, left_valley = seg(Lp1, Lp2), seg(Lv1, Lv2)
right_valley, right_peak = seg(Rv1, Rv2), seg(Rp1, Rp2)
phasing_strength = ((left_peak - left_valley) + (right_peak - right_valley)) / 2.0 \
                   if np.isfinite([left_peak, left_valley, right_peak, right_valley]).all() else float("nan")

# Emit the smoothed profile and the feature row
pd.DataFrame({"pos": pos, "mean": y}).to_csv(f"{args.out_prefix}.profile.tsv", sep="\t", index=False)
pd.DataFrame([{
  "central_mean": central_mean,
  "shoulder_mean": shoulder_mean,
  "contrast": contrast,
  "center_over_shoulder": ratio,
  "amplitude": amplitude,
  "phasing_strength": phasing_strength
}]).to_csv(f"{args.out_prefix}.features.tsv", sep="\t", index=False)
