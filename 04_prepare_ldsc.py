#!/usr/bin/env python3
"""Convert tabular GWAS summary statistics to a standard LDSC input table."""
import argparse, gzip
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
for name in ["input", "output", "chr-col", "bp-col", "snp-col", "a1-col", "a2-col", "p-col", "beta-col", "se-col", "frq-col"]:
    ap.add_argument("--" + name, required=True)
ap.add_argument("--chunksize", type=int, default=500_000)
args = ap.parse_args()
cols = [args.chr_col,args.bp_col,args.snp_col,args.a1_col,args.a2_col,args.p_col,args.beta_col,args.se_col,args.frq_col]
rename = dict(zip(cols, ["CHR","BP","SNP","A1","A2","P","BETA","SE","FRQ"]))
out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
first = True
with gzip.open(out, "wt") as handle:
    for d in pd.read_csv(args.input, sep="\t", usecols=cols, chunksize=args.chunksize):
        d = d.rename(columns=rename)
        for c in ["CHR","BP","P","BETA","SE","FRQ"]: d[c] = pd.to_numeric(d[c], errors="coerce")
        d["SNP"] = d["SNP"].fillna("").astype(str)
        fallback = d.CHR.astype("Int64").astype(str) + ":" + d.BP.astype("Int64").astype(str) + ":" + d.A1.astype(str) + ":" + d.A2.astype(str)
        d.loc[d.SNP.isin(["", "NA", "nan"]), "SNP"] = fallback
        d = d[d.CHR.between(1,22) & (d.BP > 0) & (d.P > 0) & (d.P <= 1)]
        d.to_csv(handle, sep="\t", index=False, header=first); first = False
print(f"Wrote {out}")

