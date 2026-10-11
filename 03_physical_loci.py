#!/usr/bin/env python3
"""Summarise significant variants into provisional physical lead regions.

This script is not LD clumping and must not be described as independent loci.
"""
import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
for name in ["input", "trait", "chr-col", "bp-col", "p-col", "snp-col", "a1-col", "a2-col", "beta-col", "se-col", "frq-col", "output"]:
    ap.add_argument("--" + name, required=True)
ap.add_argument("--distance", type=int, default=500_000); ap.add_argument("--chunksize", type=int, default=500_000)
args = ap.parse_args()
cols = [args.chr_col,args.bp_col,args.p_col,args.snp_col,args.a1_col,args.a2_col,args.beta_col,args.se_col,args.frq_col]
rename = dict(zip(cols, ["CHR","BP","P","SNP","A1","A2","BETA","SE","FRQ"]))
hits = []
for d in pd.read_csv(args.input, sep="\t", usecols=cols, chunksize=args.chunksize):
    d = d.rename(columns=rename)
    for c in ["CHR","BP","P","BETA","SE","FRQ"]: d[c] = pd.to_numeric(d[c], errors="coerce")
    hits.append(d[d.CHR.between(1,22) & (d.BP > 0) & (d.P < 5e-8)])
hits = pd.concat(hits, ignore_index=True).sort_values("P") if hits else pd.DataFrame(columns=rename.values())
leads = []
for _, row in hits.iterrows():
    if all(not (row.CHR == x.CHR and abs(row.BP - x.BP) <= args.distance) for x in leads): leads.append(row)
out = pd.DataFrame(leads)
if len(out):
    out.insert(0, "TRAIT", args.trait); out["LOCUS_RULE"] = f"lowest-P significant variant within plusminus{args.distance // 1000}kb; provisional physical locus"
Path(args.output).parent.mkdir(parents=True, exist_ok=True); out.to_csv(args.output, sep="\t", index=False)
print(f"Wrote {len(out)} provisional lead regions to {args.output}")

