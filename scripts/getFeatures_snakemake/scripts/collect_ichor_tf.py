#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script function: Parse ichorCNA *.params.txt files, extract (sample, tumor_fraction),
#        and write a two-column TSV: "<sample>\t<tumor_fraction>" (missing TF left blank).
# Usage:
#   python collect_ichor_tf.py --root /path/to/dir --out tf.tsv [--glob "**/*.params.txt"]
# Notes:
#   * Sample name is taken in this order: parsed from params > parent folder name > filename stem.
#   * Supports two common layouts in *.params.txt:
#       (1) A top 3-column table header "Sample  Tumor Fraction  Ploidy"
#       (2) Key-value line "Tumor Fraction: <value>" as a fallback
#   * If multiple files produce the same sample name, the last one wins (dedup behavior).
#   * Output is sorted by a natural key so e.g., ubc2 < ubc10.

import re
import argparse
from pathlib import Path

# ---------- Parse a single .params.txt ----------
def parse_params_txt(path: Path, fallback_sample: str):
    """
    Return (sample, tumor_fraction). If any is missing, it returns (sample, None).
    Sample priority: parsed from .params > parent folder name > filename stem.
    """
    txt = path.read_text(encoding="utf-8", errors="ignore")
    lines = [ln.rstrip("\n") for ln in txt.splitlines()]

    sample = None
    tf = None

    # (1) Try the top 3-column table: "Sample  Tumor Fraction  Ploidy"
    for i, ln in enumerate(lines):
        if re.search(r"^Sample\s+Tumor Fraction\s+Ploidy\s*$", ln.strip(), flags=re.IGNORECASE):
            j = i + 1
            # Skip blank lines after header
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            if j < len(lines):
                m = re.match(r"^\s*(\S+)\s+([0-9.]+)\s+([0-9.]+)\s*$", lines[j].strip())
                if m:
                    sample = m.group(1)
                    try:
                        tf = float(m.group(2))
                    except:
                        tf = None
            break

    # (2) Fallback: key-value pairs, e.g., "Tumor Fraction: 0.27"
    kv_tf = None
    lone_sample = None
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        mt = re.match(r"^Tumor Fraction:\s*([0-9.]+)\s*$", s, flags=re.IGNORECASE)
        if mt:
            try:
                kv_tf = float(mt.group(1))
            except:
                kv_tf = None
            continue
        # A single token line that could be a sample name
        if re.fullmatch(r"[A-Za-z0-9._-]+", s):
            lone_sample = lone_sample or s

    # Use KV TF if the table TF was not found
    if tf is None and kv_tf is not None:
        tf = kv_tf
    # Resolve sample name fallback chain
    if not sample:
        sample = fallback_sample or lone_sample or path.stem

    return sample, tf

# ---------- Natural sort helper (ubc1 < ubc2 < ... < ubc10) ----------
_num_pat = re.compile(r'(\d+)')
def natural_key(s: str):
    # Split into digit and non-digit parts; cast digits to int for natural ordering
    return [int(t) if t.isdigit() else t.lower() for t in _num_pat.split(s)]

# ---------- Main routine ----------
def main(root, glob_pat, out_path):
    root_path = Path(root)
    files = sorted(root_path.glob(glob_pat))
    if not files:
        print(f"[WARN] No files matched: {root_path}\\{glob_pat}")

    records = []
    for f in files:
        sample_folder = f.parent.name
        sample, tf = parse_params_txt(f, fallback_sample=sample_folder)
        records.append((str(sample), tf))

    # Deduplicate: if duplicate sample names appear, keep the last occurrence
    dedup = {}
    for name, val in records:
        dedup[name] = val

    # Sort by natural key for human-friendly ordering
    pairs = sorted(dedup.items(), key=lambda x: natural_key(x[0]))

    # Write TSV with two columns: sample \t TF (blank if TF is None)
    out_p = Path(out_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    with out_p.open("w", encoding="utf-8") as w:
        for name, val in pairs:
            w.write(f"{name}\t{'' if val is None else val}\n")

    print(f"[OK] Wrote TSV → {out_path}")
    print(f"Samples written: {len(pairs)}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Collect ichorCNA tumor fractions from *.params.txt into a two-column TSV: sample<TAB>TF.")
    ap.add_argument("--root", required=True, help="Root directory containing sample subfolders (each with a *.params.txt).")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    ap.add_argument("--glob", default="**/*.params.txt", help="Glob pattern to find parameter files (default: **/*.params.txt).")
    args = ap.parse_args()
    main(args.root, args.glob, args.out)
