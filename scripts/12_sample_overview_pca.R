#!/usr/bin/env Rscript
# Sample-level overview — PCA of the 15 pseudobulk samples, computed twice: once on
# standard VST-transformed counts (cell cycle left in), once on the cell-cycle-
# regressed residual matrix already saved by 06 (caf2_pseudobulk_cc_residuals.rds).
# Side by side, this shows directly whether removing cell cycle sharpens
# treatment-based clustering — an independent visual check on why the
# cell-cycle-independent analysis matters, not just a generic overview plot.
#
# Usage: Rscript 12_sample_overview_pca.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(DESeq2)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")

results_dir <- file.path(proj_dir, "results", CELL_LINE)
out_dir  <- file.path(results_dir, "cross_treatment_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CONDITION_LEVELS <- c("Control", "IL", "TGF", "PDGF_B", "PDGF_C")
N_TOP_VAR_GENES <- 500

# ---------------------------------------------------------------------------
# Panel A: standard analysis — VST on raw pseudobulk counts (cell cycle included)
# ---------------------------------------------------------------------------
merged <- readRDS(file.path(results_dir, "merged_filtered_ccscored.rds"))
pb <- AggregateExpression(merged, assays = "RNA", group.by = "sample", return.seurat = TRUE)
counts_mat <- as.matrix(LayerData(pb, assay = "RNA", layer = "counts"))
colnames(counts_mat) <- gsub("-", "_", colnames(counts_mat))  # undo AggregateExpression's dash sanitization

gene_filter <- read.csv(file.path(results_dir, "pseudobulk_gene_filter.csv"))
keep_genes <- gene_filter$gene[gene_filter$kept]
counts_mat <- counts_mat[keep_genes, ]

sample_meta <- merged@meta.data %>% distinct(sample, condition, replicate)
coldata <- sample_meta[match(colnames(counts_mat), sample_meta$sample), ]
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = CONDITION_LEVELS)
coldata$replicate <- factor(coldata$replicate, levels = c("r1", "r2", "r3"))
stopifnot(identical(rownames(coldata), colnames(counts_mat)))

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~condition)
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

pca_from_matrix <- function(mat, coldata, n_top = N_TOP_VAR_GENES) {
  rv <- apply(mat, 1, var)
  top <- order(rv, decreasing = TRUE)[seq_len(min(n_top, length(rv)))]
  pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
  pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df <- as.data.frame(pca$x[, 1:2])
  df$sample <- rownames(df)
  df <- cbind(df, coldata[df$sample, c("condition", "replicate")])
  list(df = df, pct_var = pct_var)
}

pca_full <- pca_from_matrix(vst_mat, coldata)

# ---------------------------------------------------------------------------
# Panel B: cell-cycle-regressed residuals (already gene-filtered to same universe)
# ---------------------------------------------------------------------------
resid_mat <- readRDS(file.path(results_dir, "pseudobulk_cc_residuals.rds"))
coldata_resid <- coldata[colnames(resid_mat), ]
pca_cc <- pca_from_matrix(resid_mat, coldata_resid)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
# Condition colors intentionally NOT set manually — using ggplot2's default
# discrete hue palette so this matches Seurat's DimPlot(group.by="condition")
# exactly (clustering_filtered/umap_overview.png), consistent across all plots.
make_pca_plot <- function(res, title) {
  df <- res$df
  hulls <- df %>% group_by(condition) %>% dplyr::slice(chull(PC1, PC2))
  ggplot(df, aes(x = PC1, y = PC2, color = condition)) +
    geom_polygon(data = hulls, aes(fill = condition), alpha = 0.15, color = NA) +
    geom_point(aes(shape = replicate), size = 3) +
    guides(fill = "none") +
    theme_minimal(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs(title = title,
         x = paste0("PC1 (", res$pct_var[1], "%)"),
         y = paste0("PC2 (", res$pct_var[2], "%)"))
}

p1 <- make_pca_plot(pca_full, "Standard analysis (VST counts, cell cycle included)")
p2 <- make_pca_plot(pca_cc, "Cell-cycle-regressed residuals")

combined <- (p1 | p2) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Sample-level overview: does removing cell cycle sharpen treatment clustering?",
    subtitle = "PCA of the 15 pseudobulk samples (top 500 most variable genes), triangles connect the 3 replicates per condition.",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) & theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "right")

out_file <- file.path(out_dir, "sample_overview_pca_before_after_cc.png")
ggsave(out_file, combined, width = 13, height = 6, dpi = 150, bg = "white")
cat("Saved:", out_file, "\n")
