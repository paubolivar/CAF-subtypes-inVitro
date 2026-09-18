#!/usr/bin/env Rscript
# Summary FGSEA barplots — one figure per treatment, one facet panel per MSigDB
# collection, top 20 pathways (10 strongest positive NES + 10 strongest negative
# NES). Bars are colored by where the pathway is significant:
#   blue         = sig_full_only          (full analysis only, cell cycle included)
#   light red    = sig_cc_independent_only (cell-cycle-independent analysis only)
#   purple       = sig_both_concordant     (robust: significant + concordant in both)
# Built from the pathway_concordance_<collection>_Control_vs_<TRT>.csv files already
# produced by 07_compare_full_vs_cc_independent.R — no FGSEA re-run needed.
#
# Usage: Rscript 08_fgsea_barplots_by_treatment.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")
comp_dir <- file.path(proj_dir, "results", CELL_LINE, "comparison_full_vs_cc_independent")

source(file.path(proj_dir, "scripts", "00_pipeline_config.R"))
source(file.path(proj_dir, "scripts", "utils_de_pathway.R"))

for (trt in TREATMENTS) {
  all_top <- list()
  for (coll in COLLECTIONS) {
    fn <- file.path(comp_dir, paste0("pathway_concordance_", coll, "_Control_vs_", trt, ".csv"))
    df <- read.csv(fn) %>%
      filter(category %in% CATEGORIES) %>%
      mutate(
        display_NES = ifelse(category == "sig_cc_independent_only", NES_cc_independent, NES_full),
        pathway_clean = clean_pathway_name(pathway, coll),
        collection = coll
      )
    top_pos <- df %>% filter(display_NES > 0) %>% slice_max(display_NES, n = 10, with_ties = FALSE)
    top_neg <- df %>% filter(display_NES < 0) %>% slice_min(display_NES, n = 10, with_ties = FALSE)
    all_top[[coll]] <- bind_rows(top_pos, top_neg)
  }
  plot_df <- bind_rows(all_top) %>%
    mutate(category = factor(category, levels = CATEGORIES)) %>%
    group_by(collection) %>%
    arrange(display_NES, .by_group = TRUE) %>%
    mutate(unique_label = factor(paste0(collection, "___", row_number()), levels = paste0(collection, "___", row_number()))) %>%
    ungroup()

  label_lookup <- setNames(plot_df$pathway_clean, as.character(plot_df$unique_label))

  p <- ggplot(plot_df, aes(x = display_NES, y = unique_label, fill = category)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey30") +
    facet_wrap(~collection, scales = "free_y", ncol = 2) +
    scale_y_discrete(labels = label_lookup) +
    scale_fill_manual(values = cat_colors, labels = cat_labels, drop = FALSE, name = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 7),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "grey90", color = NA)
    ) +
    labs(
      title = paste0("Top FGSEA pathways — Control vs ", trt),
      subtitle = "Top 10 positive + 10 negative NES per collection, colored by cell-cycle robustness",
      x = "Normalized Enrichment Score (NES)", y = NULL
    )

  out_file <- file.path(comp_dir, paste0("fgsea_barplot_", trt, ".png"))
  ggsave(out_file, p, width = 18, height = 16, dpi = 150, limitsize = FALSE, bg = "white")
  cat("Saved:", out_file, "\n")
}

# ---------------------------------------------------------------------------
# Curated companion (2026-09-14, redesigned 2026-09-14) — same padj-derived
# FDR gate as above, but pathways selected by biological theme
# (phenotype_themes / select_curated_pathways() in scripts/utils_de_pathway.R,
# applied identically across treatments) instead of an unbiased
# top-10-pos/top-10-neg by NES. Single-panel layout (theme labels in a left
# column, not faceted) via save_curated_pathway_figure() in
# scripts/utils_de_pathway.R. Additional plots only ("_curated" suffix) — the
# loop above, its output files, and their filenames are untouched.
# ---------------------------------------------------------------------------
for (trt in TREATMENTS) {
  trt_df <- bind_rows(lapply(COLLECTIONS, function(coll) {
    fn <- file.path(comp_dir, paste0("pathway_concordance_", coll, "_Control_vs_", trt, ".csv"))
    read.csv(fn) %>% mutate(collection = coll)
  }))
  curated <- select_curated_pathways(trt_df)
  if (nrow(curated) == 0) {
    cat(sprintf("No themed pathways matched for %s — skipping curated plot.\n", trt))
    next
  }

  out_file_curated <- file.path(comp_dir, paste0("fgsea_barplot_", trt, "_curated.png"))
  save_curated_pathway_figure(
    curated,
    title = paste0("Curated (CAF-phenotype) FGSEA pathways — Control vs ", trt),
    out_file = out_file_curated
  )
  cat("Saved:", out_file_curated, "\n")
}

cat("\nDone.\n")
