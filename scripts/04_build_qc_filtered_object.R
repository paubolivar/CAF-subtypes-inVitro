#!/usr/bin/env Rscript
# Build the QC-filtered, cell-cycle-scored object shared by the treatment-response
# analyses (05/06/07). Applies the SAME fixed 5-MAD nCount outlier + percent.mt>5
# filters for every cell line (not re-derived per cell line — these are fixed
# analysis-design constants, per ANALYSIS_DECISIONS.md). No automated classifier
# (e.g. scDblFinder) is run here, per standing project convention not to apply
# consequential filtering tools without first discussing them.
#
# Usage: Rscript 04_build_qc_filtered_object.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")

results_dir <- file.path(proj_dir, "results", CELL_LINE)
out_dir     <- file.path(results_dir, "qc_filtered")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged <- readRDS(file.path(results_dir, "merged_unfiltered.rds"))
cat("Loaded unfiltered merged object:", ncol(merged), "cells x", nrow(merged), "genes\n")

# --- Recompute the same 5-MAD flags as 01_qc_explore.R, self-contained ---
meta <- merged@meta.data %>%
  mutate(cell_id = rownames(.)) %>%
  group_by(sample) %>%
  mutate(
    nCount_med = median(nCount_RNA), nCount_mad = mad(nCount_RNA),
    nCount_high_outlier = nCount_RNA > (nCount_med + 5 * nCount_mad),
    nCount_low_outlier  = nCount_RNA < (nCount_med - 5 * nCount_mad),
    high_pct_mt = percent.mt > 5
  ) %>%
  ungroup()

meta$keep <- !(meta$nCount_high_outlier | meta$nCount_low_outlier | meta$high_pct_mt)

# --- Audit table: before/after per sample ---
audit <- meta %>%
  group_by(sample, condition, replicate) %>%
  summarise(
    n_cells_before             = n(),
    n_removed_high_ncount      = sum(nCount_high_outlier),
    n_removed_low_ncount       = sum(nCount_low_outlier),
    n_removed_high_mt          = sum(high_pct_mt & !nCount_high_outlier & !nCount_low_outlier),
    n_cells_after              = sum(keep),
    pct_removed                = round(100 * (1 - n_cells_after / n_cells_before), 2),
    .groups = "drop"
  ) %>%
  arrange(condition, replicate)

write.csv(audit, file.path(out_dir, "qc_filter_audit.csv"), row.names = FALSE)
cat("\n--- QC filter audit (per sample) ---\n")
print(as.data.frame(audit), row.names = FALSE)
cat("\nTotal cells before:", nrow(meta), " after:", sum(meta$keep),
    " (", round(100 * (1 - sum(meta$keep) / nrow(meta)), 2), "% removed)\n", sep = "")

# --- Diagnostic plots: before vs after ---
meta$qc_status <- ifelse(meta$keep, "kept", "removed")

p1 <- ggplot(meta, aes(x = sample, y = nCount_RNA, fill = qc_status)) +
  geom_violin(scale = "width") + scale_y_log10() +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "nCount_RNA per sample — kept vs removed", x = NULL, y = "UMI count")

p2 <- ggplot(meta, aes(x = sample, y = percent.mt, fill = qc_status)) +
  geom_violin(scale = "width") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "percent.mt per sample — kept vs removed", x = NULL, y = "% mitochondrial")

ggsave(file.path(out_dir, "qc_filter_before_after_violins.png"), (p1 / p2), width = 12, height = 8, dpi = 150)

cells_after <- meta %>% filter(keep) %>% count(sample, condition, replicate)
p3 <- ggplot(cells_after, aes(x = sample, y = n, fill = condition)) +
  geom_col() + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "Cells per sample after QC filtering", y = "n cells")
ggsave(file.path(out_dir, "qc_filter_cells_per_sample.png"), p3, width = 8, height = 5, dpi = 150)

# --- Subset to kept cells ---
keep_cells <- meta$cell_id[meta$keep]
merged <- subset(merged, cells = keep_cells)
cat("\nAfter QC filtering:", ncol(merged), "cells x", nrow(merged), "genes\n")

# --- Join layers, normalize, cell-cycle score (same call pattern as 02_cluster_explore.R) ---
merged <- JoinLayers(merged)
merged <- NormalizeData(merged, verbose = FALSE)

s_genes   <- intersect(cc.genes.updated.2019$s.genes,   rownames(merged))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(merged))
cat("\nCell-cycle gene overlap with panel: S genes", length(s_genes), "/", length(cc.genes.updated.2019$s.genes),
    " | G2M genes", length(g2m_genes), "/", length(cc.genes.updated.2019$g2m.genes), "\n")
merged <- CellCycleScoring(merged, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

cat("\nPhase distribution after filtering:\n")
print(table(merged$Phase))

# No ScaleData/PCA/UMAP/clustering here — not needed for pseudobulk DGE, and avoids
# implying a clustering decision that hasn't been made (clusters 15/16 still undecided).

out_rds <- file.path(results_dir, "merged_filtered_ccscored.rds")
saveRDS(merged, out_rds)
cat("\nDone. Saved:", out_rds, "\n")
cat("Outputs in:", out_dir, "\n")
