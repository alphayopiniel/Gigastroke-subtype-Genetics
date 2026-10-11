#!/usr/bin/env bash
# Run non-functional SuSiE fine-mapping with a UK Biobank British-ancestry LD block.
# Requires PolyFun from its upstream repository. Do not commit the downloaded LD block.
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 POLYFUN_DIR SUMSTATS CHR START END SAMPLE_SIZE" >&2
  exit 1
fi

polyfun_dir=$1
sumstats=$2
chr=$3
start=$4
end=$5
n=$6

# Set these paths outside this repository, for example in a local shell profile.
: "${UKBB_LD_DIR:?Set UKBB_LD_DIR to the directory containing UKBB LD npz/gz blocks}"
: "${OUTFILE:?Set OUTFILE to a results TSV path outside the repository}"

ld_prefix="${UKBB_LD_DIR}/chr${chr}_${start}_${end}"
python "${polyfun_dir}/finemapper.py" \
  --sumstats "${sumstats}" \
  --chr "${chr}" \
  --start "${start}" \
  --end "${end}" \
  --n "${n}" \
  --ld "${ld_prefix}" \
  --method susie \
  --non-funct \
  --max-num-causal 10 \
  --susie-resvar 1 \
  --memory 8 \
  --out "${OUTFILE}"

