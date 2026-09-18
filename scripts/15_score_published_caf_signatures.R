#!/usr/bin/env Rscript
# Score published CAF gene signatures against the QC-filtered CAF2 object via
# Seurat::AddModuleScore, after harmonizing mouse gene symbols to human via
# Human_mouse_orthology_v1.txt (orthology_confidence==1 only, per-gene: direct
# match to our panel first, else mouse->human ortholog lookup, else drop).
#
# Produces one 3-panel figure per signature:
#   1. Violin + boxplot (quantiles) of per-cell score by condition
#   2. Overlaid density of per-cell score by condition
#   3. UMAP feature plot, faceted by condition (so a small high-scoring
#      subpopulation within one treatment is visible, not hidden in one combined blob)
#
# Runs across all 72 kept signatures (Cords/Elyada/Gao/Mechta/Pietras/Wu/KPlab_curated — see
# script 13 for which sources/signatures are included; Ruth/Ostman/Ohlund excluded
# per user decision). Design was fine-tuned on Wu alone first before this full run.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")

results_dir <- file.path(proj_dir, "results", CELL_LINE)
sig_dir  <- file.path(proj_dir, "data", "published_data", "CAF_gene_signatures")
out_dir  <- file.path(results_dir, "caf_signature_scoring")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SOURCES_TO_RUN <- c("Cords", "Elyada", "Gao", "Mechta", "Pietras", "Wu", "KPlab_curated")  # all kept sources (Ruth/Ostman/Ohlund excluded per user decision)
CONDITION_LEVELS <- c("Control", "IL", "TGF", "PDGF_B", "PDGF_C")

set.seed(1)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
merged <- readRDS(file.path(results_dir, "merged_filtered_ccscored_umap.rds"))
merged$condition <- factor(merged$condition, levels = CONDITION_LEVELS)
cat("Loaded:", ncol(merged), "cells x", nrow(merged), "genes\n")

signatures <- read.csv(file.path(sig_dir, "CAF_gene_signatures_combined.csv")) %>%
  filter(source %in% SOURCES_TO_RUN)
cat("Signatures to run:", paste(unique(signatures$signature), collapse = ", "), "\n")

# --- Orthology table (mouse -> human), confidence==1 only ---
ortho <- read.delim(file.path(proj_dir, "data", "published_data", "Human_mouse_orthology_v1.txt"),
                     sep = "\t", stringsAsFactors = FALSE)
ortho <- ortho[ortho$orthology_confidence == 1, ]
mouse_to_human <- setNames(ortho$Homsa_gene_name, ortho$Musmu_gene_name)
cat("Orthology table (confidence==1):", nrow(ortho), "mouse->human mappings\n")

# ---------------------------------------------------------------------------
# 2. Harmonize gene symbols per signature (per-gene: direct match -> ortholog -> drop)
# ---------------------------------------------------------------------------
panel_genes <- rownames(merged)

harmonize_signature <- function(genes) {
  direct <- genes %in% panel_genes
  mapped <- ifelse(direct, genes, mouse_to_human[genes])
  mapped_valid <- !is.na(mapped) & mapped %in% panel_genes
  data.frame(
    original_gene = genes,
    final_gene = ifelse(mapped_valid, mapped, NA_character_),
    matched_direct = direct,
    matched_ortholog = !direct & mapped_valid,
    stringsAsFactors = FALSE
  )
}

log_rows <- list()
sig_gene_lists <- list()
for (sig in unique(signatures$signature)) {
  genes <- signatures %>% filter(signature == sig) %>% arrange(rank) %>% pull(gene)
  h <- harmonize_signature(genes)
  final_genes <- h$final_gene[!is.na(h$final_gene)]
  sig_gene_lists[[sig]] <- final_genes
  log_rows[[sig]] <- data.frame(
    signature = sig,
    n_genes_total = length(genes),
    n_matched_direct = sum(h$matched_direct),
    n_matched_ortholog = sum(h$matched_ortholog),
    n_dropped = sum(is.na(h$final_gene)),
    n_used = length(final_genes)
  )
}
harmonization_log <- bind_rows(log_rows)
write.csv(harmonization_log, file.path(out_dir, "harmonization_log.csv"), row.names = FALSE)
cat("\n--- Harmonization log ---\n")
print(harmonization_log, row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Score + plot per signature
# ---------------------------------------------------------------------------
TREATMENTS <- setdiff(CONDITION_LEVELS, "Control")
summary_rows <- list()

for (sig in names(sig_gene_lists)) {
  genes <- sig_gene_lists[[sig]]
  if (length(genes) < 2) {
    cat("Skipping", sig, "- fewer than 2 genes available after harmonization\n")
    next
  }

  merged <- AddModuleScore(merged, features = list(genes), name = "sig_score_tmp", ctrl = 100)
  score_values <- merged@meta.data[["sig_score_tmp1"]]
  merged@meta.data[["sig_score_tmp1"]] <- NULL

  df <- data.frame(
    score = score_values,
    condition = merged$condition,
    UMAP_1 = Embeddings(merged, "umap")[, 1],
    UMAP_2 = Embeddings(merged, "umap")[, 2]
  )
  df_umap <- df[order(df$score), ]  # plot high scorers last/on top, not buried under the majority

  # --- Per-treatment summary vs Control (Wilcoxon rank-sum on per-cell scores) ---
  control_scores <- df$score[df$condition == "Control"]
  mean_control <- mean(control_scores)
  for (trt in TREATMENTS) {
    trt_scores <- df$score[df$condition == trt]
    wt <- wilcox.test(trt_scores, control_scores)
    summary_rows[[paste(sig, trt)]] <- data.frame(
      signature = sig,
      source = sub("_.*", "", sig),
      treatment = trt,
      mean_control = mean_control,
      mean_treatment = mean(trt_scores),
      delta = mean(trt_scores) - mean_control,
      wilcox_pvalue = wt$p.value
    )
  }

  n_used <- length(genes)
  n_total <- harmonization_log$n_genes_total[harmonization_log$signature == sig]
  subtitle_txt <- paste(strwrap(paste0(
    n_used, "/", n_total, " genes used after harmonization: ", paste(genes, collapse = ", ")
  ), width = 140), collapse = "\n")

  # Condition colors intentionally NOT set manually here — using ggplot2's default
  # discrete hue palette so this matches Seurat's DimPlot(group.by="condition")
  # exactly (results/clustering/umap_overview.png), consistent across all plots.
  p_violin <- ggplot(df, aes(x = condition, y = score, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.6, linewidth = 0.4) +
    guides(fill = "none") +
    theme_minimal(base_size = 9) +
    labs(title = "Score by condition", x = NULL, y = "Module score") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  p_density <- ggplot(df, aes(x = score, color = condition, fill = condition)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    theme_minimal(base_size = 9) +
    labs(title = "Score density by condition", x = "Module score", y = "Density") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.position = "right")

  # Diverging scale centered at this signature's own MEDIAN score, not a fixed 0 —
  # AddModuleScore output isn't reliably centered at 0 in practice (e.g. myCAFs
  # skews positive here), so a fixed-0 midpoint can push almost all cells onto one
  # side of the scale and wash out contrast. Centering on the median guarantees the
  # color scale always splits the data into below-/above-typical regardless of
  # where the raw distribution sits. Points ordered ascending by score above so
  # high-scoring cells draw on top, not buried under the majority.
  score_midpoint <- median(df_umap$score)
  p_umap <- ggplot(df_umap, aes(x = UMAP_1, y = UMAP_2, color = score)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = score_midpoint, name = "score") +
    facet_wrap(~condition, nrow = 1) +
    theme_minimal(base_size = 9) +
    labs(title = "UMAP feature plot, by condition") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          strip.background = element_rect(fill = "grey90", color = NA))

  combined <- (p_violin | p_density) / p_umap +
    plot_annotation(
      title = sig,
      subtitle = subtitle_txt,
      theme = theme(plot.background = element_rect(fill = "white", color = NA))
    ) & theme(plot.background = element_rect(fill = "white", color = NA))

  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", sig)
  out_file <- file.path(out_dir, paste0("signature_score_", safe_name, ".png"))
  ggsave(out_file, combined, width = 14, height = 9, dpi = 150, bg = "white")
  cat("Saved:", out_file, "\n")
}

# ---------------------------------------------------------------------------
# 4. Cross-signature summary table (feeds 17_cross_signature_bubble.R)
# ---------------------------------------------------------------------------
signature_treatment_summary <- bind_rows(summary_rows)
signature_treatment_summary$wilcox_padj <- p.adjust(signature_treatment_summary$wilcox_pvalue, method = "BH")
summary_out <- file.path(out_dir, "signature_treatment_summary.csv")
write.csv(signature_treatment_summary, summary_out, row.names = FALSE)
cat("Saved:", summary_out, "\n")

cat("\nDone.\n")
