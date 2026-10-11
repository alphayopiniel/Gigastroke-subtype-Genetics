#!/usr/bin/env python3
"""Position/allele harmonisation for two regional summary-statistic files.

Input files must already have CHR, BP, SNP, A1, A2, BETA, SE and P columns.
The output retains only exact matches or simple strand-compatible allele swaps.
"""
import argparse
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--left", required=True); ap.add_argument("--right", required=True)
ap.add_argument("--chr", required=True, type=int); ap.add_argument("--start", required=True, type=int); ap.add_argument("--end", required=True, type=int)
ap.add_argument("--output", required=True)
args = ap.parse_args()
need = ["CHR","BP","SNP","A1","A2","BETA","SE","P"]
def read(path, suffix):
    d = pd.read_csv(path, sep="\t", usecols=need)
    d = d[(d.CHR == args.chr) & d.BP.between(args.start, args.end)].copy()
    return d.rename(columns={c: f"{c}_{suffix}" for c in need if c not in ["CHR", "BP"]})

left, right = read(args.left, "LEFT"), read(args.right, "RIGHT")
d = left.merge(right, on=["CHR","BP"], how="inner")
same = (d.A1_LEFT == d.A1_RIGHT) & (d.A2_LEFT == d.A2_RIGHT)
swap = (d.A1_LEFT == d.A2_RIGHT) & (d.A2_LEFT == d.A1_RIGHT)
d = d[same | swap].copy()
d["RIGHT_BETA_ALIGNED_TO_LEFT_A1"] = d.BETA_RIGHT.where(same.loc[d.index], -d.BETA_RIGHT)
d["ALLELE_ACTION"] = "same"; d.loc[swap.loc[d.index], "ALLELE_ACTION"] = "swapped"
d = d.sort_values(["CHR","BP"])
Path(args.output).parent.mkdir(parents=True, exist_ok=True); d.to_csv(args.output, sep="\t", index=False)
print(f"Aligned {len(d)} variants; wrote {args.output}")

