#!/usr/bin/env Rscript
# Part A — FULL treatment-response analysis (cell cycle INCLUDED, not adjusted for).
# Pseudobulk DGE (Control vs each of 4 treatments) via DESeq2 on raw summed counts,
# plus FGSEA pathway analysis across 5 MSigDB collections. Biological replicate
# (sample) is the unit of inference (ANALYSIS_DECISIONS.md, 2026-09-09).
#
# Usage: Rscript 05_pseudobulk_treatment_response_full.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(DESeq2)
  library(edgeR)
  library(msigdbr)
  library(fgsea)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")

results_dir <- file.path(proj_dir, "results", CELL_LINE)
out_dir     <- file.path(results_dir, "pseudobulk_full")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(proj_dir, "scripts", "00_pipeline_config.R"))
source(file.path(proj_dir, "scripts", "utils_de_pathway.R"))
set.seed(1)

merged <- readRDS(file.path(results_dir, "merged_filtered_ccscored.rds"))
cat("Loaded QC-filtered, CC-scored object:", ncol(merged), "cells x", nrow(merged), "genes\n")

# ---------------------------------------------------------------------------
# 1. Pseudobulk aggregation (raw counts, summed per sample)
# ---------------------------------------------------------------------------
pb <- AggregateExpression(merged, assays = "RNA", group.by = "sample", return.seurat = TRUE)
counts_mat <- as.matrix(LayerData(pb, assay = "RNA", layer = "counts"))
# AggregateExpression() sanitizes group.by identities (underscores -> dashes) in the
# resulting column names; restore the original "sample" values (which only ever
# contain underscores as separators, never dashes) so they match sample metadata.
colnames(counts_mat) <- gsub("-", "_", colnames(counts_mat))
cat("\nPseudobulk matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")
stopifnot(all(counts_mat >= 0), all(counts_mat == round(counts_mat)))
cat("Sanity check passed: pseudobulk counts are non-negative integers.\n")

# --- coldata, one row per sample, matching counts_mat column order ---
sample_meta <- merged@meta.data %>%
  distinct(sample, condition, replicate)
coldata <- sample_meta[match(colnames(counts_mat), sample_meta$sample), ]
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = CONDITION_LEVELS)
coldata$replicate <- factor(coldata$replicate, levels = REPLICATE_LEVELS)
stopifnot(identical(rownames(coldata), colnames(counts_mat)))
cat("\n--- coldata ---\n")
print(coldata)

# ---------------------------------------------------------------------------
# 2. Gene filtering (shared with Part B for a fair comparison)
# ---------------------------------------------------------------------------
keep <- filterByExpr(counts_mat, group = coldata$condition)
cat("\nGenes kept by filterByExpr:", sum(keep), "/", length(keep), "\n")
gene_filter_df <- data.frame(gene = rownames(counts_mat), kept = keep)
write.csv(gene_filter_df, file.path(results_dir, "pseudobulk_gene_filter.csv"), row.names = FALSE)
counts_mat <- counts_mat[keep, ]

# ---------------------------------------------------------------------------
# 3. DESeq2
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~condition)
dds <- DESeq(dds)
cat("\nresultsNames(dds):\n")
print(resultsNames(dds))

all_res <- list()
for (trt in TREATMENTS) {
  coef_name <- paste0("condition_", trt, "_vs_Control")
  res <- results(dds, contrast = c("condition", trt, "Control"))
  if (coef_name %in% resultsNames(dds)) {
    res_shrunk <- tryCatch(
      lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
      error = function(e) {
        cat("apeglm failed for", coef_name, "- falling back to ashr:", conditionMessage(e), "\n")
        lfcShrink(dds, contrast = c("condition", trt, "Control"), type = "ashr", quiet = TRUE)
      }
    )
  } else {
    cat("Coefficient", coef_name, "not found in resultsNames(dds); using ashr on contrast instead.\n")
    res_shrunk <- lfcShrink(dds, contrast = c("condition", trt, "Control"), type = "ashr", quiet = TRUE)
  }
  # stat/pvalue/padj come from the unshrunk `res`; log2FoldChange/lfcSE from shrunk
  df <- as.data.frame(res_shrunk)
  df$stat   <- res$stat
  df$pvalue <- res$pvalue
  df$padj   <- res$padj
  df$gene <- rownames(df)
  df$contrast <- paste0("Control_vs_", trt)
  df <- df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "contrast")]
  df <- df[order(df$padj), ]
  write.csv(df, file.path(out_dir, paste0("deseq2_Control_vs_", trt, ".csv")), row.names = FALSE)
  all_res[[trt]] <- df
  n_sig <- sum(!is.na(df$padj) & df$padj < SIG_ALPHA & abs(df$log2FoldChange) >= fc_cutoff("deseq2"))
  cat(sprintf("\n%s: %d genes tested, %d significant at padj<%s\n", trt, nrow(df), n_sig, SIG_ALPHA))
}
deseq2_all <- bind_rows(all_res)
write.csv(deseq2_all, file.path(out_dir, "deseq2_all_contrasts.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. FGSEA — 5 MSigDB collections (shared runner: scripts/utils_de_pathway.R)
# ---------------------------------------------------------------------------
fgsea_all_df <- run_fgsea_all(all_res, TREATMENTS, stat_col = "stat", out_dir = out_dir)

# ---------------------------------------------------------------------------
# 5. Plots
# ---------------------------------------------------------------------------
volc_plots <- lapply(TREATMENTS, function(trt) {
  df <- all_res[[trt]]
  df$sig <- !is.na(df$padj) & df$padj < SIG_ALPHA & abs(df$log2FoldChange) >= fc_cutoff("deseq2")
  build_volcano_plot(df, x_col = "log2FoldChange", y_pval_col = "pvalue", padj_col = "padj",
                      title = paste0("Control vs ", trt), xlab = "log2FC")
})
ggsave(file.path(out_dir, "volcano_all_contrasts.png"),
       wrap_plots(volc_plots, ncol = 2), width = 10, height = 9, dpi = 150)

write_fgsea_dotplots(fgsea_all_df, out_dir, title_suffix = "full analysis, cell cycle included")

# ---------------------------------------------------------------------------
# 6. Brief interpretation
# ---------------------------------------------------------------------------
lines <- c(
  "# Full treatment-response analysis (cell cycle included)",
  "",
  "This is the FULL analysis (`pseudobulk_full/`), referred to as `full` in",
  "the comparison outputs. Pseudobulk DGE (DESeq2, raw summed counts per biological",
  "replicate) comparing each treatment to Control. Cell-cycle composition",
  "differences between samples are NOT adjusted for here — see",
  "`pseudobulk_cc_independent/` (the `cc_independent` analysis) and",
  "`comparison_full_vs_cc_independent/` for the cell-cycle-independent view",
  "and the comparison between the two.",
  ""
)
for (trt in TREATMENTS) {
  df <- all_res[[trt]]
  cutoff <- fc_cutoff("deseq2")
  n_sig <- sum(!is.na(df$padj) & df$padj < SIG_ALPHA & abs(df$log2FoldChange) >= cutoff)
  top_up <- df %>% filter(!is.na(padj), padj < SIG_ALPHA, log2FoldChange >= cutoff) %>% slice_max(log2FoldChange, n = 20)
  top_dn <- df %>% filter(!is.na(padj), padj < SIG_ALPHA, log2FoldChange <= -cutoff) %>% slice_min(log2FoldChange, n = 20)
  lines <- c(lines,
    paste0("## Control vs ", trt),
    paste0("- ", n_sig, " genes meeting padj<", SIG_ALPHA, " and |log2FC|>=", cutoff, "."),
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
writeLines(lines, file.path(out_dir, "interpretation_full.md"))

cat("\nDone. Outputs in:", out_dir, "\n")
