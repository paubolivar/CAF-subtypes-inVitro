#!/usr/bin/env Rscript
# Cross-treatment DGE heatmap — per treatment, the union of the top DESeq2/full
# DEGs and the top limma/cc-independent DEGs (each method ranked and selected
# independently, using the SAME established DEG criteria the rest of the
# pipeline uses), pooled across all 4 treatments into one gene list, shown
# across all 4 treatments at once. Two side-by-side panels (full analysis /
# cc-independent analysis) sharing the same clustered gene order, built
# directly from the gene concordance CSV already produced by
# 07_compare_full_vs_cc_independent.R — no re-analysis needed.
#
# GENE SELECTION (changed 2026-09-14 — union of both methods, not
# cc-independent alone): the previous version selected purely by top
# |logFC_cc_independent| among genes significant in the cc-independent
# analysis, with DESeq2/full shown alongside for context only but never
# itself driving selection. That made the heatmap asymmetric between the two
# methods. It's now built from BOTH methods independently:
#   1. top TOP_N_PER_TREATMENT DESeq2 DEGs per treatment (gene_all$sig_full,
#      i.e. padj_full < SIG_ALPHA AND |log2FoldChange_full| >=
#      fc_cutoff("deseq2") — 07's precomputed flag), ranked by
#      abs(log2FoldChange_full);
#   2. top TOP_N_PER_TREATMENT limma DEGs per treatment
#      (gene_all$sig_cc_independent, i.e. adj.P.Val_cc_independent <
#      SIG_ALPHA AND |logFC_cc_independent| >= fc_cutoff("limma")), ranked by
#      abs(logFC_cc_independent);
#   3. the UNION (not intersection) of those two lists, per treatment — a
#      gene selected by only one method is intentionally kept, since showing
#      where the two methods agree AND where they disagree is the point of
#      this figure;
#   4. the per-treatment unions are pooled (union + dedup) across all 4
#      treatments into the final plotted gene set.
# Both panels show this identical gene set in identical (clustered) order —
# see the empirical verification printed below and reported back after each
# run. A gene shown in a panel is NOT required to be significant in that
# panel's own analysis; the `*` annotation (using the same sig_full/
# sig_cc_independent columns as selection, not a separately-invented
# padj-only rule) marks exactly which cells individually meet the DEG
# criteria.
#
# SCALE CAVEAT (shown on the plot too): the full analysis's log2FoldChange (DESeq2,
# base-2 log fold change of counts) and the cc-independent analysis's logFC (limma,
# natural-log-scale residuals after cell-cycle regression) are NOT on the same
# numeric scale — they agree in sign much of the time but not necessarily in
# magnitude. Each panel therefore gets its own independently-scaled color legend; do
# not visually compare color intensity between the two panels.
#
# Usage: Rscript 10_cross_treatment_dge_heatmap.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
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
comp_dir <- file.path(results_dir, "comparison_full_vs_cc_independent")
out_dir  <- file.path(results_dir, "cross_treatment_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(proj_dir, "scripts", "00_pipeline_config.R"))

TOP_N_PER_TREATMENT <- 15

gene_all <- read.csv(file.path(comp_dir, "gene_concordance_all_contrasts.csv"))

# ---------------------------------------------------------------------------
# Per-treatment: top-N DESeq2 DEGs UNION top-N limma DEGs, each independently
# ranked and filtered on the established (padj + effect-size) DEG criteria
# already computed by 07 (`sig_full`, `sig_cc_independent`) — not a
# locally-invented threshold.
# ---------------------------------------------------------------------------
deseq2_top <- gene_all %>%
  filter(sig_full) %>%
  group_by(contrast) %>%
  slice_max(abs(log2FoldChange_full), n = TOP_N_PER_TREATMENT, with_ties = FALSE) %>%
  ungroup()

limma_top <- gene_all %>%
  filter(sig_cc_independent) %>%
  group_by(contrast) %>%
  slice_max(abs(logFC_cc_independent), n = TOP_N_PER_TREATMENT, with_ties = FALSE) %>%
  ungroup()

per_trt_report <- lapply(TREATMENTS, function(trt) {
  ctr <- paste0("Control_vs_", trt)
  d_genes <- deseq2_top %>% filter(contrast == ctr) %>% pull(gene)
  l_genes <- limma_top  %>% filter(contrast == ctr) %>% pull(gene)
  data.frame(
    treatment       = trt,
    deseq2_top15    = length(d_genes),
    limma_top15     = length(l_genes),
    common          = length(intersect(d_genes, l_genes)),
    deseq2_only     = length(setdiff(d_genes, l_genes)),
    limma_only      = length(setdiff(l_genes, d_genes)),
    treatment_union = length(union(d_genes, l_genes))
  )
}) %>% bind_rows()
cat("\n=== Per-treatment DESeq2 vs limma top-", TOP_N_PER_TREATMENT, " overlap ===\n", sep = "")
print(per_trt_report, row.names = FALSE)

top_genes <- unique(c(deseq2_top$gene, limma_top$gene))
n_deseq2_any <- length(intersect(top_genes, unique(deseq2_top$gene)))
n_limma_any  <- length(intersect(top_genes, unique(limma_top$gene)))
n_both_any   <- length(intersect(unique(deseq2_top$gene), unique(limma_top$gene)))
cat("\nTotal unique pooled genes (final heatmap gene set):", length(top_genes), "\n")
cat("Selected via DESeq2 in >=1 treatment:", n_deseq2_any, "\n")
cat("Selected via limma in >=1 treatment:", n_limma_any, "\n")
cat("Selected via BOTH methods in >=1 treatment (not necessarily the same treatment):", n_both_any, "\n")

plot_data <- gene_all %>% filter(gene %in% top_genes)
stopifnot(nrow(plot_data) == length(top_genes) * length(TREATMENTS))  # every pooled gene has a real row for every treatment
stopifnot(!anyNA(plot_data$log2FoldChange_full), !anyNA(plot_data$logFC_cc_independent))

# --- Build wide matrices for clustering ---
mat_full <- plot_data %>% select(gene, contrast, log2FoldChange_full) %>%
  pivot_wider(names_from = contrast, values_from = log2FoldChange_full) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()
mat_cc <- plot_data %>% select(gene, contrast, logFC_cc_independent) %>%
  pivot_wider(names_from = contrast, values_from = logFC_cc_independent) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()
mat_cc <- mat_cc[rownames(mat_full), ]  # ensure same row order pre-clustering

# --- Cluster rows on the cc-independent (trustworthy) values ---
hc <- hclust(dist(mat_cc), method = "average")
gene_order <- rownames(mat_cc)[hc$order]

contrast_labels <- setNames(paste0("Control_vs_", TREATMENTS), paste0("Control_vs_", TREATMENTS))

plot_long <- plot_data %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    contrast = factor(contrast, levels = paste0("Control_vs_", TREATMENTS), labels = TREATMENTS)
    # sig_full / sig_cc_independent are used AS-IS from gene_all (07's
    # established padj+effect-size DEG criteria) for the `*` annotation below
    # — no local, looser redefinition.
  )

cat("\nBoth panels use the same `gene` factor (identical genes, identical order):",
    identical(levels(plot_long$gene), gene_order), "\n")
cat("Every plotted value traces to gene_concordance_all_contrasts.csv's own",
    "log2FoldChange_full / logFC_cc_independent columns (no averaging, no substitution) — ",
    "verified structurally above (row count check, NA check).\n")

base_theme <- theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  )

p_full <- ggplot(plot_long, aes(x = contrast, y = gene, fill = log2FoldChange_full)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(sig_full, "*", "")), color = "black", size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "log2FC\n(full)") +
  labs(title = "Full analysis", x = NULL, y = NULL) +
  base_theme

p_cc <- ggplot(plot_long, aes(x = contrast, y = gene, fill = logFC_cc_independent)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(sig_cc_independent, "*", "")), color = "black", size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "logFC\n(cc-independent,\nnatural-log residuals)") +
  labs(title = "Cell-cycle-independent analysis", x = NULL, y = NULL) +
  base_theme + theme(axis.text.y = element_blank())

combined <- (p_full | p_cc) +
  plot_annotation(
    title = "Top DESeq2 (full) + limma (cc-independent) DEGs across all 4 treatments",
    subtitle = paste(strwrap(paste0(
                       "Per treatment: union of the top ", TOP_N_PER_TREATMENT, " DESeq2 DEGs (by |log2FoldChange_full|) and top ", TOP_N_PER_TREATMENT,
                       " limma DEGs (by |logFC_cc_independent|), each independently significant per the established DEG criteria; pooled (", length(gene_order),
                       " unique genes) and clustered by cc-independent pattern. ",
                       "* = meets the established DEG criteria (padj + effect-size) in that specific analysis for that treatment — a gene can be shown without an asterisk if it was selected via the other method or another treatment. ",
                       "Color scales differ between panels (see script 10 comments) — compare direction/significance, not raw color intensity, across panels."
                     ), width = 150), collapse = "\n"),
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.background = element_rect(fill = "white", color = NA))

out_file <- file.path(out_dir, "cross_treatment_dge_heatmap.png")
ggsave(out_file, combined, width = 12, height = max(8, length(gene_order) * 0.22), dpi = 150, limitsize = FALSE, bg = "white")
cat("Saved:", out_file, "\n")
