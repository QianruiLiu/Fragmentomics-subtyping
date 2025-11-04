#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Aggregate Griffin site-level coverage metrics into a wide sample×feature table.
#        Scans per-sample directories, prefers rep2 when multiple reps exist, expands each
#        site into columns (AD_/NE_ prefixes) and exports a TSV. Supports raw and GC-corrected inputs.
# Usage:
#   python collect_griffin_sitewise.py -b /path/to/griffin/results -o merged/griffin_sitewise_features.tsv --metrics mean_coverage,central_coverage,amplitude [--only_gccorr] [--keep_dir_name] [--no_recal_suffix]
# Notes:
#   * Input files are expected as *.coverage.tsv or *.GC_corrected.coverage.tsv inside each sample dir.
#   * Site names must contain "AD-ATAC" or "NE-ATAC" to route into AD_/NE_ columns.
#   * Column names emitted as: AD_<site_name>_<metric>_(GCcorr|raw).

import os, re, argparse
from glob import glob
from collections import defaultdict, Counter
from typing import Optional, Dict, List
import pandas as pd

# Directories to ignore while scanning sample folders
PLOT_DIRS = {"plots", "plot", "logs", "tmp", "temp", "figs", "images"}

def is_ad(site): return "AD-ATAC" in str(site)
def is_ne(site): return "NE-ATAC" in str(site)

def parse_rep_from_name(name: str) -> Optional[int]:
    """Extract replication number from a directory/file name, e.g., '_rep2_' -> 2; return None if absent."""
    m = re.search(r'[_\-\.]rep(\d+)(?:[_\-.]|$)', name, flags=re.I)
    return int(m.group(1)) if m else None

def base_stem_dir(dir_name: str) -> str:
    """
    A "base stem" to group duplicate dirs across reps by stripping '_rep\\d+'.
    Example: '173_1_rep2_LuCaP' and '173_1_rep1_LuCaP' -> '173_1_LuCaP'
    """
    return re.sub(r'([_\-])?rep\d+(?=[_\-]|$)', '', dir_name, flags=re.I)

def normalize_sample_name_from_dir(dir_name: str, keep_original: bool, rep_for_name: Optional[int], append_recal: bool) -> str:
    """Sample naming policy: keep original directory name by default; optionally standardize to LuCaP_XX_repN_recal."""
    if keep_original:
        return dir_name
    # Standardized style: LuCaP_XX_repN_recal (best-effort digit extraction)
    m = re.search(r'(\d+)', dir_name)
    base = f"LuCaP_{m.group(1)}" if m else dir_name
    if rep_for_name is not None:
        base += f"_rep{rep_for_name}"
    if append_recal and not base.endswith("_recal"):
        base += "_recal"
    return base

def pick_metric_columns(df: pd.DataFrame):
    """
    Identify three metric columns (no averaging here; used to locate value columns):
      mean_coverage    <- Window-Mean / window_mean / WindowMean / mean_coverage / ...
      central_coverage <- Central-Mean / central_mean / CentralMean / ...
      amplitude        <- any name containing 'amplit'
    """
    cols = list(df.columns)
    norm_map = {c.lower().replace("-", "_"): c for c in cols}

    def find_first(keys):
        for k in keys:
            key = k.lower().replace("-", "_")
            if key in norm_map:
                return norm_map[key]
        return None

    mean_cov = find_first(["mean_coverage", "window_mean", "windowmean", "window-mean"])
    central_cov = find_first(["central_coverage", "central_mean", "centralmean", "central-mean"])
    amp_candidates = [c for c in cols if re.search(r"amplit", c, re.I)]
    amp_col = None
    if amp_candidates:
        exact = [c for c in amp_candidates if re.fullmatch(r"(?i)amplitude", c)]
        amp_col = exact[0] if exact else amp_candidates[0]

    return mean_cov, central_cov, amp_col

def read_one_file_sitewise(path: str, suffix: str, take_metrics: List[str]) -> Dict[str, float]:
    """
    Read a single *.coverage.tsv and expand by site_name into columns (no averaging).
    Column naming:  AD_<site>_<metric>_<suffix>   or  NE_<site>_<metric>_<suffix>
      e.g.,         AD_AR_TFBS_000123_mean_coverage_GCcorr
    take_metrics may include: mean_coverage / central_coverage / amplitude (any subset).
    """
    df = pd.read_csv(path, sep="\t")
    # Align site_name column
    if "site_name" not in df.columns:
        alt = [c for c in df.columns if c.lower() == "site_name"]
        if alt:
            df = df.rename(columns={alt[0]: "site_name"})
        else:
            return {}

    mean_cov, central_cov, amp_col = pick_metric_columns(df)

    # Keep only requested metrics that are present
    metric_map = {
        "mean_coverage": mean_cov,
        "central_coverage": central_cov,
        "amplitude": amp_col,
    }
    metric_cols = {m: c for m, c in metric_map.items() if m in take_metrics and c is not None}
    if not metric_cols:
        return {}

    # Subset to necessary columns
    need = ["site_name"] + list(metric_cols.values())
    df = df[need].copy()

    out = {}

    # Expand into AD_/NE_ columns without averaging
    for tag, mask_fn in (("AD", is_ad), ("NE", is_ne)):
        sub = df[df["site_name"].apply(mask_fn)].copy()
        if sub.empty:
            continue
        # Sanitize site_name so it is safe as a column label
        sub["site_name"] = (
            sub["site_name"].astype(str)
            .str.replace(r"\s+", "_", regex=True)
            .str.replace(r"[^A-Za-z0-9_.\-]", "_", regex=True)
        )
        for mkey, col in metric_cols.items():
            # One column per site
            for _, row in sub.iterrows():
                val = pd.to_numeric(row[col], errors="coerce")
                if pd.notna(val):
                    colname = f"{tag}_{row['site_name']}_{mkey}_{suffix}"
                    out[colname] = float(val)

    return out

def main():
    ap = argparse.ArgumentParser(
        description="Aggregate Griffin coverage (site-level) into a wide sample×feature table; prefer rep2 per sample stem."
    )
    ap.add_argument(
        "-b", "--base_dir",
        default=r"\\wsl$\Ubuntu\home\student\1_Griffin\snakemakes\griffin_nucleosome_profiling\results",
        help="Root directory containing per-sample subfolders."
    )
    ap.add_argument(
        "-o", "--output",
        default="merged/griffin_sitewise_features.tsv",
        help="Output TSV path (rows=samples; columns=site×metric×suffix)."
    )
    ap.add_argument(
        "--metrics",
        default="mean_coverage,central_coverage,amplitude",
        help="Comma-separated metrics to collect. Choices: mean_coverage,central_coverage,amplitude"
    )
    ap.add_argument(
        "--only_gccorr",
        action="store_true",
        help="Use only *.GC_corrected.coverage.tsv (by default, collect both raw and GCcorr)."
    )
    ap.add_argument(
        "--keep_dir_name",
        action="store_true",
        help="Keep original directory name as sample ID (default). If not set, standardize to LuCaP_XX_repN_recal."
    )
    ap.add_argument(
        "--no_recal_suffix",
        action="store_true",
        help="When standardizing name, do not append '_recal'. Only effective if --keep_dir_name is not set."
    )
    args = ap.parse_args()

    take_metrics = [x.strip() for x in args.metrics.split(",") if x.strip()]
    for m in take_metrics:
        if m not in {"mean_coverage", "central_coverage", "amplitude"}:
            raise SystemExit(f"[ERR] Unsupported metric: {m}")

    base = args.base_dir
    if not os.path.exists(base):
        raise SystemExit("Base dir not found: {}".format(base))

    # 1) Enumerate sample directories; prefer rep2 at the directory level
    all_dirs = [d for d in glob(os.path.join(base, "*")) if os.path.isdir(d)]
    all_dirs = [d for d in all_dirs if os.path.basename(d).lower() not in PLOT_DIRS]
    if not all_dirs:
        raise SystemExit("No sample directories under: {}".format(base))

    grouped: Dict[str, List[str]] = defaultdict(list)
    for sdir in all_dirs:
        dn = os.path.basename(sdir)
        stem = base_stem_dir(dn)
        grouped[stem].append(sdir)

    selected_dirs: List[str] = []
    for stem, members in grouped.items():
        # If a rep2 dir exists, select rep2; otherwise pick the only one,
        # or choose the highest rep number if present; else lexicographically last.
        rep2 = [d for d in members if parse_rep_from_name(os.path.basename(d)) == 2]
        if rep2:
            chosen = sorted(rep2)[-1]
        else:
            if len(members) == 1:
                chosen = members[0]
            else:
                with_rep = [(parse_rep_from_name(os.path.basename(d)) or -1, d) for d in members]
                chosen = sorted(with_rep)[-1][1]
        selected_dirs.append(chosen)

    if not selected_dirs:
        raise SystemExit("No directories selected after rep2 filtering.")

    # 2) For each selected directory, collect TSVs and expand site-wise
    sample_to_feats = {}
    for sdir in sorted(selected_dirs):
        raw_dir = os.path.basename(sdir)

        gc_paths = [p for p in glob(os.path.join(sdir, "*.GC_corrected.coverage.tsv"))
                    if re.search(r"\.GC_corrected\.coverage\.tsv$", p)]
        raw_paths = [] if args.only_gccorr else [
            p for p in glob(os.path.join(sdir, "*.coverage.tsv"))
            if not re.search(r"\.GC_corrected\.coverage\.tsv$", p)
        ]
        all_paths = gc_paths + raw_paths
        if not all_paths:
            print(f"[INFO] No coverage TSV in {raw_dir}")
            continue

        # Use replication number from directory name when constructing standardized sample name
        rep_for_name = parse_rep_from_name(raw_dir)

        sample_name = normalize_sample_name_from_dir(
            raw_dir,
            keep_original=args.keep_dir_name or True,  # default True: keep original dir name
            rep_for_name=rep_for_name,
            append_recal=(not args.no_recal_suffix)
        )

        bag = {}
        for pth in sorted(all_paths):
            suffix = "GCcorr" if pth.endswith(".GC_corrected.coverage.tsv") else "raw"
            kv = read_one_file_sitewise(pth, suffix, take_metrics)
            # Later entries override earlier ones for same column (usually suffix differs so no collision)
            bag.update(kv)

        if not bag:
            print(f"[INFO] No usable site-wise metrics in {raw_dir}")
            continue

        sample_to_feats[sample_name] = bag

    if not sample_to_feats:
        raise SystemExit("No data aggregated. Check paths, filenames, and metrics.")

    # 3) Assemble DataFrame: rows=samples; columns=site×metric×suffix
    all_cols = set()
    for d in sample_to_feats.values():
        all_cols.update(d.keys())

    # Column ordering: AD first → NE second; metrics mean→central→amplitude; GCcorr first → raw second; then by site name
    def col_key(c: str):
        tag = 0 if c.startswith("AD_") else 1
        metric_rank = 3
        if "_mean_coverage_" in c:
            metric_rank = 0
        elif "_central_coverage_" in c:
            metric_rank = 1
        elif "_amplitude_" in c:
            metric_rank = 2
        suffix_rank = 0 if c.endswith("_GCcorr") else 1
        site_part = re.sub(r'^(AD|NE)_', '', c)
        site_part = re.sub(r'_(mean_coverage|central_coverage|amplitude)_(GCcorr|raw)$', '', site_part)
        return (tag, metric_rank, suffix_rank, site_part, c)

    cols_sorted = sorted(all_cols, key=col_key)

    df = pd.DataFrame.from_dict(sample_to_feats, orient="index")
    df = df.reindex(columns=cols_sorted)
    df = df.sort_index()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, sep="\t", index=True, header=True)
    print(f"[OK] wrote {args.output} with shape {df.shape}")
    print("Rows = samples (directory name or standardized name), Columns = AD/NE_<site>_<metric>_(GCcorr|raw)")
    print("Dedup policy: keep only rep2 directory for the same base stem; if no rep2, pick best by heuristic.")
if __name__ == "__main__":
    main()
