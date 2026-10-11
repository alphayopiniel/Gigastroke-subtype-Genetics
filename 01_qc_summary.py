#!/usr/bin/env python3
"""Streaming QC summary for a tab-separated GWAS summary-statistics file."""
import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True)
ap.add_argument("--trait", required=True)
ap.add_argument("--p-col", required=True)
ap.add_argument("--output", required=True)
ap.add_argument("--chunksize", type=int, default=500_000)
args = ap.parse_args()

total = valid = invalid = gws = 0
for chunk in pd.read_csv(args.input, sep="\t", usecols=[args.p_col], chunksize=args.chunksize):
    values = pd.to_numeric(chunk[args.p_col], errors="coerce")
    ok = (values > 0) & (values <= 1)
    total += len(values); valid += int(ok.sum()); invalid += int((~ok).sum())
    gws += int((values[ok] < 5e-8).sum())

text = "\n".join([
    f"Trait: {args.trait}", f"Total variants: {total}",
    f"Valid P-values: {valid}", f"Invalid or missing P-values: {invalid}",
    f"Genome-wide significant variants (P < 5e-8): {gws}", ""
])
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
Path(args.output).write_text(text)
print(text, end="")

