#!/usr/bin/env Rscript
# Curated (investigator-selected) FGSEA barplots — a companion to the unbiased
# top-20-by-NES plots from 08_fgsea_barplots_by_treatment.R, NOT a replacement.
#
# Rationale: the unbiased top-N plots include statistically significant but
# biologically uninformative-for-this-question pathways (e.g. generic antiviral
# defense gene sets, unrelated tissue gene sets like cardiac/neuronal/sensory terms
# that happen to reach significance via shared housekeeping genes). This script
# instead selects, from the SAME already-computed, already-multiple-testing-corrected
# FGSEA results, only pathways matching an explicit relevance keyword rule tied to
# the project's actual biological question (CAF identity/heterogeneity axes:
# inflammatory/immune signaling, ECM/collagen, myofibroblast/contractility,
# proliferation). This is a display/curation choice, not a new statistical test — no
# p-values are recomputed, and the unbiased plot remains the primary record. The
# exact include/exclude keyword rule used is printed on the plot itself and in the
# console so the selection is auditable.
#
# Usage: Rscript 09_fgsea_barplot_selected.R [CELL_LINE]   (default: CAF2)
# The relevance rules below are about treatment biology (IL/TGF/PDGF), not cell-line
# biology, so they apply unchanged to any cell line run through this same design.

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

MAX_PER_COLLECTION <- 15  # cap so panels stay readable even if many terms match

# `selection_rules` (per-treatment biological relevance rule, auditable, edit
# in scripts/utils_de_pathway.R to extend to other treatments) is defined once
# there and shared with 11_cross_treatment_pathway_bubble.R so the two plots'
# "curated" pathway sets cannot silently desync.
TREATMENTS <- names(selection_rules)

for (trt in TREATMENTS) {
  rule <- selection_rules[[trt]]
  inc_pat <- paste(rule$include, collapse = "|")
  exc_pat <- paste(rule$exclude, collapse = "|")
  cat(sprintf("\n[%s] include if name matches: %s\n", trt, inc_pat))
  cat(sprintf("[%s] excluded even if matched: %s\n", trt, exc_pat))

  all_sel <- list()
  for (coll in COLLECTIONS) {
    fn <- file.path(comp_dir, paste0("pathway_concordance_", coll, "_Control_vs_", trt, ".csv"))
    df <- read.csv(fn) %>%
      filter(category %in% CATEGORIES) %>%
      filter(grepl(inc_pat, pathway, ignore.case = TRUE)) %>%
      filter(!grepl(exc_pat, pathway, ignore.case = TRUE)) %>%
      mutate(
        display_NES = ifelse(category == "sig_cc_independent_only", NES_cc_independent, NES_full),
        pathway_clean = clean_pathway_name(pathway, coll),
        collection = coll
      ) %>%
      slice_max(abs(display_NES), n = MAX_PER_COLLECTION, with_ties = FALSE)
    cat(sprintf("  %s: %d pathways selected\n", coll, nrow(df)))
    all_sel[[coll]] <- df
  }
  plot_df <- bind_rows(all_sel)

  if (nrow(plot_df) == 0) {
    cat(sprintf("No pathways matched the relevance rule for %s — skipping plot.\n", trt))
    next
  }

  plot_df <- plot_df %>%
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
      strip.background = element_rect(fill = "grey90", color = NA),
      plot.caption = element_text(hjust = 0, size = 7, color = "grey30")
    ) +
    labs(
      title = paste0("Selected FGSEA pathways — Control vs ", trt, " (investigator-curated)"),
      subtitle = "Significant pathways selected for relevance to this treatment's expected biology; NOT an unbiased top-N (see fgsea_barplot_*.png for that)",
      x = "Normalized Enrichment Score (NES)", y = NULL,
      caption = paste(strwrap(paste0(
        "Selection rule (see script 09 for exact regex) — included: ", rule$include_readable,
        ". Excluded even if significant: ", rule$exclude_readable, ". ",
        "All pathways shown remain individually significant (padj<0.05) and multiple-testing-corrected; only the DISPLAY subset was curated, not the statistics."
      ), width = 130), collapse = "\n")
    )

  out_file <- file.path(comp_dir, paste0("fgsea_barplot_", trt, "_selected.png"))
  ggsave(out_file, p, width = 16, height = 12, dpi = 150, limitsize = FALSE, bg = "white")
  cat("Saved:", out_file, "\n")
}

cat("\nDone.\n")
