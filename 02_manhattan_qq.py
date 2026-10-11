#!/usr/bin/env python3
"""Create Manhattan and QQ plots; labels are display-only, not LD clumps."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True); ap.add_argument("--trait", required=True)
ap.add_argument("--chr-col", required=True); ap.add_argument("--bp-col", required=True)
ap.add_argument("--p-col", required=True); ap.add_argument("--snp-col", required=True)
ap.add_argument("--output-dir", required=True); ap.add_argument("--top-n", type=int, default=10)
ap.add_argument("--label-distance", type=int, default=2_000_000)
ap.add_argument("--chunksize", type=int, default=500_000)
args = ap.parse_args(); out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

lengths = [248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468]
offsets, ticks, cumulative = {}, [], 0
for chrom, length in enumerate(lengths, 1):
    offsets[chrom] = cumulative; ticks.append(cumulative + length / 2); cumulative += length

use = [args.chr_col, args.bp_col, args.p_col, args.snp_col]
def clean(d):
    d = d.rename(columns={args.chr_col:"CHR", args.bp_col:"BP", args.p_col:"P", args.snp_col:"SNP"})
    for c in ["CHR", "BP", "P"]: d[c] = pd.to_numeric(d[c], errors="coerce")
    return d[d.CHR.between(1, 22) & (d.BP > 0) & (d.P > 0) & (d.P <= 1)].copy()

pvals, candidates = [], []
for d in pd.read_csv(args.input, sep="\t", usecols=use, chunksize=args.chunksize):
    d = clean(d); pvals.append(d.P.to_numpy()); candidates.append(d.nsmallest(max(100, args.top_n * 10), "P"))
labels = []
for _, row in pd.concat(candidates).sort_values("P").iterrows():
    if all(not (int(row.CHR) == int(x.CHR) and abs(row.BP - x.BP) <= args.label_distance) for x in labels): labels.append(row)
    if len(labels) >= args.top_n: break

plt.figure(figsize=(18, 7))
for d in pd.read_csv(args.input, sep="\t", usecols=use, chunksize=args.chunksize):
    d = clean(d)
    for chrom in range(1, 23):
        x = d[d.CHR == chrom]
        if len(x): plt.scatter(x.BP + offsets[chrom], -np.log10(x.P), s=.15, c="#4C78A8" if chrom % 2 else "#F58518", linewidths=0)
plt.axhline(-np.log10(5e-8), color="red", linestyle="--", linewidth=1)
for i, row in enumerate(labels):
    name = str(row.SNP) if pd.notna(row.SNP) else f"chr{int(row.CHR)}:{int(row.BP)}"
    plt.annotate(name, (row.BP + offsets[int(row.CHR)], -np.log10(row.P)), xytext=(0, 8 + 8 * (i % 3)), textcoords="offset points", ha="center", fontsize=8, rotation=45)
plt.xticks(ticks, range(1, 23)); plt.xlabel("Chromosome"); plt.ylabel("-log10(P)"); plt.title(f"Manhattan plot: {args.trait}")
plt.tight_layout(); plt.savefig(out / f"{args.trait}_manhattan.png", dpi=300); plt.close()

p = np.sort(np.concatenate(pvals)); n = len(p); take = np.unique(np.linspace(0, n - 1, min(n, 500_000), dtype=int))
expected, observed = -np.log10((take + 1) / (n + 1)), -np.log10(p[take]); lim = 1.03 * max(expected.max(), observed.max())
plt.figure(figsize=(7, 7)); plt.scatter(expected, observed, s=.3, alpha=.6, color="#4C78A8"); plt.plot([0, lim], [0, lim], "--", color="red")
plt.xlim(0, lim); plt.ylim(0, lim); plt.xlabel("Expected -log10(P)"); plt.ylabel("Observed -log10(P)"); plt.title(f"QQ plot: {args.trait}")
plt.tight_layout(); plt.savefig(out / f"{args.trait}_qq.png", dpi=300)

