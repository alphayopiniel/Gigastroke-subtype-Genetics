#!/usr/bin/env python3
"""Combine CES/LAS provisional leads into a preregistered replication manifest."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--ces", required=True); ap.add_argument("--las", required=True); ap.add_argument("--output", required=True)
args = ap.parse_args()
ces, las = pd.read_csv(args.ces, sep="\t"), pd.read_csv(args.las, sep="\t")
for label, d in [("CES", ces), ("LAS", las)]:
    d["TRAIT"] = label
    d["DISCOVERY_OR"] = np.exp(pd.to_numeric(d["BETA"], errors="coerce"))
    d["EXPECTED_EFFECT_DIRECTION"] = np.where(d["BETA"] > 0, "A1 increases subtype risk", "A1 decreases subtype risk")
    d["REPLICATION_REQUIREMENT"] = "Non-overlapping CES/LAS cohort; same A1 allele and effect direction; two-sided P < 0.005"
    d["DISCOVERY_SOURCE"] = "GIGASTROKE standard subtype GWAS"
out = pd.concat([ces, las], ignore_index=True)
core = ["TRAIT","CHR","BP","SNP","A1","A2","BETA","SE","P","FRQ","DISCOVERY_OR","EXPECTED_EFFECT_DIRECTION","LOCUS_RULE","DISCOVERY_SOURCE","REPLICATION_REQUIREMENT"]
out = out[[c for c in core if c in out.columns]]
Path(args.output).parent.mkdir(parents=True, exist_ok=True); out.to_csv(args.output, sep="\t", index=False)
print(f"Wrote {len(out)} replication candidates to {args.output}")

