# Subtype Specific Genetic Architecture of Ischaemic Stroke

Reproducible workflow accompanying **"Subtype-Specific Genetic Architecture of Ischaemic Stroke Reveals Distinct Cardioembolic, Atherosclerotic, and Small-Vessel Patterns."**

## Scope

The primary analyses use standard GIGASTROKE summary statistics for cardioembolic stroke (CES), large artery stroke (LAS), and small vessel stroke (SVS). The workflow produces:

1. QC summaries and Manhattan/QQ plots;
2. provisional physical lead-locus tables;
3. LDSC-ready files and commands for subtype genetic correlations;
4. harmonised regional inputs for CES-AF and LAS-CAD analyses;
5. replication manifests for the seven CES and three LAS provisional lead variants.

The repository contains **no** GIGASTROKE, UK Biobank LD, 1000 Genomes LD, or replication data. Obtain and use every input according to its own terms. Do not upload controlled-access data, individual-level data, LD matrices, or derived files that disclose restricted data.

## Important scientific interpretation

* Physical lead regions are a transparent reporting summary, not LD clumping or fine-mapping.
* A high single-signal colocalisation posterior does not establish a shared causal variant when multi-signal diagnostics are discordant or non-evaluable.
* The CES-AF chromosome 4 result was assessed using a 503-sample 1000 Genomes EUR LD reference; dense-LD SuSiE-coloc supported distinct regional signals.
* The LAS-CAD chromosome 7 multi-signal analysis was non-evaluable because LAS did not form a high-purity credible set with either 1000 Genomes EUR or UK Biobank British-ancestry LD.
* The provisional lead variants require replication in a non-overlapping cohort with CES/LAS phenotyping.

## Requirements

* Python 3.10 or later
* `numpy`, `pandas`, `matplotlib`
* LDSC and the appropriate European LD-score reference files for genetic-correlation analyses
* R with `data.table`, `coloc`, and `susieR` only for the optional regional colocalisation script

```bash
python -m pip install -r requirements.txt
```

## Expected input layout

```text
data/
  CES_sumstats.tsv.gz
  LAS_sumstats.tsv.gz
  SVS_sumstats.tsv.gz
  AF_sumstats.tsv.gz
  CAD_sumstats.tsv.gz
reference/                         # ignored by Git
results/                           # ignored by Git
```

The scripts expect tab-separated, headered summary statistics. Supply the actual column names using command-line options; no assumption is made that all sources use the same schema.

## Workflow

Set paths once:

```bash
DATA=data
RESULTS=results
LDSC=/path/to/ldsc
HM3=/path/to/eur_w_ld_chr/w_hm3.snplist
REF=/path/to/eur_w_ld_chr
mkdir -p "$RESULTS"
```

### 1. QC and discovery plots

```bash
for trait in CES LAS SVS; do
  python scripts/01_qc_summary.py --input "$DATA/${trait}_sumstats.tsv.gz" --trait "$trait" --p-col P --output "$RESULTS/${trait}_qc.txt"
  python scripts/02_manhattan_qq.py --input "$DATA/${trait}_sumstats.tsv.gz" --trait "$trait" --chr-col CHR --bp-col BP --p-col P --snp-col SNP --output-dir "$RESULTS/plots"
done
```

### 2. Provisional physical locus summaries

```bash
python scripts/03_physical_loci.py --input "$DATA/CES_sumstats.tsv.gz" --trait CES --chr-col CHR --bp-col BP --p-col P --snp-col SNP --a1-col A1 --a2-col A2 --beta-col BETA --se-col SE --frq-col FRQ --output "$RESULTS/CES_provisional_lead_loci.tsv"
python scripts/03_physical_loci.py --input "$DATA/LAS_sumstats.tsv.gz" --trait LAS --chr-col CHR --bp-col BP --p-col P --snp-col SNP --a1-col A1 --a2-col A2 --beta-col BETA --se-col SE --frq-col FRQ --output "$RESULTS/LAS_provisional_lead_loci.tsv"
```

### 3. LDSC preparation and genetic correlations

```bash
for trait in CES LAS SVS AF CAD; do
  python scripts/04_prepare_ldsc.py --input "$DATA/${trait}_sumstats.tsv.gz" --output "$RESULTS/${trait}_for_ldsc.tsv.gz" --chr-col CHR --bp-col BP --snp-col SNP --a1-col A1 --a2-col A2 --p-col P --beta-col BETA --se-col SE --frq-col FRQ
done

cd "$LDSC"
python munge_sumstats.py --sumstats "$RESULTS/CES_for_ldsc.tsv.gz" --snp SNP --a1 A1 --a2 A2 --p P --frq FRQ --signed-sumstats BETA,0 --N-cas 12790 --N-con 1503898 --merge-alleles "$HM3" --out "$RESULTS/CES_ldsc"
# Repeat munging with the correct effective sample size for every trait.
python ldsc.py --rg "$RESULTS/CES_ldsc.sumstats.gz,$RESULTS/AF_ldsc.sumstats.gz" --ref-ld-chr "$REF/" --w-ld-chr "$REF/" --out "$RESULTS/CES_AF_rg"
```

### 4. Prepare regional inputs

```bash
python scripts/05_harmonise_region.py --left "$RESULTS/CES_for_ldsc.tsv.gz" --right "$RESULTS/AF_for_ldsc.tsv.gz" --chr 4 --start 109000000 --end 112000000 --output "$RESULTS/CES_AF_chr4_harmonised.tsv"
python scripts/05_harmonise_region.py --left "$RESULTS/LAS_for_ldsc.tsv.gz" --right "$RESULTS/CAD_for_ldsc.tsv.gz" --chr 7 --start 18000001 --end 21000001 --output "$RESULTS/LAS_CAD_chr7_harmonised.tsv"
```

Use the saved harmonised output with an LD matrix matched to the chosen reference panel. `scripts/06_coloc_susie.R` is deliberately diagnostic: it refuses to claim a shared causal variant when credible-set information is absent.

### 4b. UK Biobank LD sensitivity fine-mapping for LAS-CAD

Download the required UK Biobank LD block under its access terms, install PolyFun from its upstream repository, and keep both outside this repository. Then run each trait separately:

```bash
export UKBB_LD_DIR=/path/to/UKBB_LD
export OUTFILE="$RESULTS/LAS_chr7_UKBB_susie.tsv"
bash scripts/08_ukbb_susie.sh /path/to/polyfun "$RESULTS/LAS_for_ldsc.tsv.gz" 7 18000001 21000001 1514702

export OUTFILE="$RESULTS/CAD_chr7_UKBB_susie.tsv"
bash scripts/08_ukbb_susie.sh /path/to/polyfun "$RESULTS/CAD_for_ldsc.tsv.gz" 7 18000001 21000001 1165360
```

Audit credible-set membership and PIP concentration before running, or interpreting, multi-signal colocalisation. A trait with no high-purity credible set makes the signal-pair comparison non-evaluable.

### 5. Create the replication manifest

```bash
python scripts/07_replication_manifest.py --ces "$RESULTS/CES_provisional_lead_loci.tsv" --las "$RESULTS/LAS_provisional_lead_loci.tsv" --output "$RESULTS/replication_manifest_CES_LAS.tsv"
```

## Citation and acknowledgements

Please cite the GIGASTROKE primary publication: Mishra A, Malik R, Hachiya T, et al. *Nature*. 2022;611:115-123. doi:10.1038/s41586-022-05165-3. Cite LDSC, coloc/SuSiE, and every trait-specific GWAS source used. Include all consortium-mandated acknowledgements.

See [code provenance](docs/PROVENANCE.md) for what is, and is not, included in this public repository.

