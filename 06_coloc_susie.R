#!/usr/bin/env Rscript
# Diagnostic regional SuSiE-coloc analysis.
# Usage: Rscript scripts/06_coloc_susie.R harmonised.tsv ld_matrix.rds left_n right_n output_prefix
suppressPackageStartupMessages({ library(data.table); library(susieR); library(coloc) })
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 5) stop("Usage: 06_coloc_susie.R harmonised.tsv ld_matrix.rds left_n right_n output_prefix")
input <- fread(argv[1]); R <- readRDS(argv[2]); n_left <- as.numeric(argv[3]); n_right <- as.numeric(argv[4]); prefix <- argv[5]
if (nrow(input) != nrow(R) || ncol(R) != nrow(R)) stop("LD matrix dimensions do not match harmonised variants")
z_left <- input$BETA_LEFT / input$SE_LEFT
z_right <- input$RIGHT_BETA_ALIGNED_TO_LEFT_A1 / input$SE_RIGHT
s_left <- susie_rss(z = z_left, R = R, n = n_left, estimate_residual_variance = FALSE)
s_right <- susie_rss(z = z_right, R = R, n = n_right, estimate_residual_variance = FALSE)
result <- coloc.susie(s_left, s_right)
saveRDS(list(susie_left=s_left, susie_right=s_right, coloc=result), paste0(prefix, ".rds"))
if (is.null(result$summary) || !nrow(result$summary)) {
  fwrite(data.table(note="non-evaluable: no comparable credible-set pair"), paste0(prefix, "_summary.tsv"), sep="\t")
} else fwrite(as.data.table(result$summary), paste0(prefix, "_summary.tsv"), sep="\t")

