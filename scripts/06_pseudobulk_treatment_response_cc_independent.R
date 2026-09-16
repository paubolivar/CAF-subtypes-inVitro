#!/usr/bin/env Rscript
# Part B — CELL-CYCLE-INDEPENDENT treatment-response analysis.
# S.Score/G2M.Score are regressed out of each cell's log-normalized expression
# (ScaleData, do.scale=FALSE, do.center=FALSE -> residuals, not z-scores), then
# residuals are mean-aggregated per sample into a pseudobulk matrix and tested with
# limma (NOT DESeq2 — these are continuous residuals, not integer counts, so a
# negative-binomial count model is not appropriate here).
#
# IMPORTANT CAVEAT (see also interpretation_cc_independent.md and Part C): because
# Part A uses DESeq2 on raw counts and Part B uses limma on continuous CC-regressed
# residuals, the two parts are not on identical statistical footing. Differences
# between A and B partly reflect this methodological difference, not purely
# cell-cycle removal.
#
# Usage: Rscript 06_pseudobulk_treatment_response_cc_independent.R [CELL_LINE]
# (default: CAF2)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(limma)
  library(msigdbr)
  library(fgsea)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")

results_dir <- file.path(proj_dir, "results", CELL_LINE)
out_dir     <- file.path(results_dir, "pseudobulk_cc_independent")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(proj_dir, "scripts", "00_pipeline_config.R"))
source(file.path(proj_dir, "scripts", "utils_de_pathway.R"))
set.seed(1)

merged <- readRDS(file.path(results_dir, "merged_filtered_ccscored.rds"))
cat("Loaded QC-filtered, CC-scored object:", ncol(merged), "cells x", nrow(merged), "genes\n")
stopifnot(all(c("S.Score", "G2M.Score") %in% colnames(merged@meta.data)))

# Restrict regression up front to the same gene universe used in Part A
# (edgeR::filterByExpr) — these are the only genes analyzed downstream anyway, and
# this is an expression-level filter, not a variable-feature/HVG selection, so it
# does not bias the pathway analysis background.
gene_filter <- read.csv(file.path(results_dir, "pseudobulk_gene_filter.csv"))
regress_genes <- intersect(gene_filter$gene[gene_filter$kept], rownames(merged))
cat("Regressing cell-cycle out of", length(regress_genes), "genes (Part A's filterByExpr set)\n")

# ---------------------------------------------------------------------------
# 1. Regress out cell-cycle scores per cell (Part A's gene set, not just HVGs)
# ---------------------------------------------------------------------------
# Every gene shares the identical design matrix (~S.Score+G2M.Score), so per-gene
# lm() fitting (Seurat's ScaleData(vars.to.regress=...) default) is unnecessary:
# the OLS residuals for ALL genes at once are a single matrix projection. This is
# mathematically identical to Seurat's default linear-model regression (residuals =
# observed - fitted from `lm(gene ~ S.Score + G2M.Score)`), computed via one BLAS
# matrix multiplication instead of ~14k sequential per-gene fits (which took >3h
# single-threaded, and ran into the exact same wall via future-parallelized
# ScaleData due to whole-object serialization to each worker).
data_slot <- LayerData(merged, assay = "RNA", layer = "data")[regress_genes, , drop = FALSE]
Y <- as.matrix(Matrix::t(data_slot))  # cells x genes, dense
X <- model.matrix(~ S.Score + G2M.Score, data = merged@meta.data)  # cells x 3 (intercept + 2 covariates)
stopifnot(identical(rownames(X), colnames(data_slot)))

beta   <- solve(crossprod(X), crossprod(X, Y))  # 3 x genes
fitted <- X %*% beta                            # cells x genes
resid  <- Y - fitted                            # cells x genes

resid_mat <- t(resid)
rownames(resid_mat) <- regress_genes
colnames(resid_mat) <- colnames(data_slot)
cat("\nResidual matrix (post CC-regression):", nrow(resid_mat), "genes x", ncol(resid_mat), "cells\n")
stopifnot(nrow(resid_mat) == length(regress_genes))
rm(Y, fitted, resid, data_slot); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 2. Mean-aggregate residuals per sample (mean, not sum — continuous, not counts)
# ---------------------------------------------------------------------------
samples <- sort(unique(merged$sample))
pb_resid <- sapply(samples, function(s) {
  Matrix::rowMeans(resid_mat[, merged$sample == s, drop = FALSE])
})
rownames(pb_resid) <- rownames(resid_mat)
cat("\nPseudobulk residual matrix:", nrow(pb_resid), "genes x", ncol(pb_resid), "samples\n")

saveRDS(pb_resid, file.path(results_dir, "pseudobulk_cc_residuals.rds"))
stopifnot(setequal(rownames(pb_resid), regress_genes))  # confirms Part A gene universe match

# --- coldata ---
sample_meta <- merged@meta.data %>% distinct(sample, condition, replicate)
coldata <- sample_meta[match(colnames(pb_resid), sample_meta$sample), ]
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = CONDITION_LEVELS)
coldata$replicate <- factor(coldata$replicate, levels = REPLICATE_LEVELS)
stopifnot(identical(rownames(coldata), colnames(pb_resid)))

# ---------------------------------------------------------------------------
# 4. limma
# ---------------------------------------------------------------------------
design <- model.matrix(~condition, data = coldata)
cat("\nDesign matrix columns:", paste(colnames(design), collapse = ", "), "\n")
fit <- lmFit(pb_resid, design)
fit <- eBayes(fit)

all_res <- list()
for (trt in TREATMENTS) {
  coef_name <- paste0("condition", trt)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$contrast <- paste0("Control_vs_", trt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "contrast")]
  tt <- tt[order(tt$adj.P.Val), ]
  write.csv(tt, file.path(out_dir, paste0("limma_Control_vs_", trt, ".csv")), row.names = FALSE)
  all_res[[trt]] <- tt
  n_sig <- sum(!is.na(tt$adj.P.Val) & tt$adj.P.Val < SIG_ALPHA & abs(tt$logFC) >= fc_cutoff("limma"), na.rm = TRUE)
  cat(sprintf("\n%s: %d genes tested, %d significant at adj.P.Val<%s\n", trt, nrow(tt), n_sig, SIG_ALPHA))
}
limma_all <- bind_rows(all_res)
write.csv(limma_all, file.path(out_dir, "limma_all_contrasts.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. FGSEA — same 5 MSigDB collections as Part A, ranked by limma's moderated t
#    (shared runner: scripts/utils_de_pathway.R)
# ---------------------------------------------------------------------------
fgsea_all_df <- run_fgsea_all(all_res, TREATMENTS, stat_col = "t", out_dir = out_dir)

# ---------------------------------------------------------------------------
# 6. Plots
# ---------------------------------------------------------------------------
volc_plots <- lapply(TREATMENTS, function(trt) {
  df <- all_res[[trt]]
  # limma logFC is on the residual expression scale here; fc_cutoff("limma")
  # converts the shared LOG2FC_CUTOFF onto that natural-log scale (see
  # 00_pipeline_config.R) so the same nominal effect-size gate applies.
  df$sig <- !is.na(df$adj.P.Val) &
            df$adj.P.Val < SIG_ALPHA &
            abs(df$logFC) >= fc_cutoff("limma")
  build_volcano_plot(df, x_col = "logFC", y_pval_col = "P.Value", padj_col = "adj.P.Val",
                      title = paste0("Control vs ", trt, " (CC-independent)"),
                      xlab = "logFC (limma residual scale)")
})
ggsave(file.path(out_dir, "volcano_all_contrasts.png"),
       wrap_plots(volc_plots, ncol = 2), width = 10, height = 9, dpi = 150)

write_fgsea_dotplots(fgsea_all_df, out_dir, title_suffix = "cell-cycle-independent analysis")

# ---------------------------------------------------------------------------
# 7. Brief interpretation
# ---------------------------------------------------------------------------
lines <- c(
  "# Cell-cycle-independent treatment-response analysis",
  "",
  "This is the CELL-CYCLE-INDEPENDENT analysis",
  "(`pseudobulk_cc_independent/`), referred to as `cc_independent` in the",
  "comparison outputs. S.Score/G2M.Score were regressed out of each cell's",
  "log-normalized expression (mathematically equivalent to Seurat ScaleData,",
  "do.scale=FALSE, do.center=FALSE), residuals were mean-aggregated per biological",
  "replicate, and tested with limma (moderated t-statistics).",
  "",
  "**Statistical-footing caveat**: the FULL analysis (`pseudobulk_full/`,",
  "`full` in the comparison outputs) used DESeq2 on raw summed pseudobulk counts",
  "(negative-binomial GLM). This cell-cycle-independent analysis uses limma on",
  "continuous CC-regressed residuals (moderated linear model on Gaussian-like data).",
  "These are different statistical frameworks with different null models,",
  "variance-shrinkage mechanisms, and power characteristics. Differences between the",
  "full and cell-cycle-independent analyses (see",
  "`comparison_full_vs_cc_independent/`) therefore partly reflect this",
  "methodological gap, not purely cell-cycle removal — treat 'CC-driven' conclusions",
  "there as suggestive, not definitive.",
  ""
)
for (trt in TREATMENTS) {
  df <- all_res[[trt]]
  cutoff <- fc_cutoff("limma")
  n_sig <- sum(!is.na(df$adj.P.Val) & df$adj.P.Val < SIG_ALPHA & abs(df$logFC) >= cutoff)
  top_up <- df %>% filter(!is.na(adj.P.Val), adj.P.Val < SIG_ALPHA, logFC >= cutoff) %>% slice_max(logFC, n = 20)
  top_dn <- df %>% filter(!is.na(adj.P.Val), adj.P.Val < SIG_ALPHA, logFC <= -cutoff) %>% slice_min(logFC, n = 20)
  lines <- c(lines,
    paste0("## Control vs ", trt),
    paste0("- ", n_sig, " genes meeting adj.P.Val<", SIG_ALPHA, " and equivalent |log2FC|>=", LOG2FC_CUTOFF, " cutoff."),
    paste0("- Top upregulated: ", paste(top_up$gene, collapse = ", ")),
    paste0("- Top downregulated: ", paste(top_dn$gene, collapse = ", ")),
    ""
  )
  for (coll_name in unique(fgsea_all_df$collection)) {
    fsub <- fgsea_all_df %>% filter(collection == coll_name, contrast == paste0("Control_vs_", trt), !is.na(padj), padj < SIG_ALPHA)
    if (nrow(fsub) == 0) next
    top_path <- fsub %>% slice_max(abs(NES), n = 3)
    lines <- c(lines, paste0("- ", coll_name, " top pathways: ", paste(top_path$pathway, collapse = "; ")))
  }
  lines <- c(lines, "")
}
writeLines(lines, file.path(out_dir, "interpretation_cc_independent.md"))

cat("\nDone. Outputs in:", out_dir, "\n")
