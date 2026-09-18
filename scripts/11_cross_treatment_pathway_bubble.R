#!/usr/bin/env Rscript
# Cross-treatment pathway bubble plot — pools each treatment's curated
# biologically-relevant pathways (same relevance rules as
# 09_fgsea_barplot_selected.R) into one combined pathway list, then shows every
# pooled pathway across ALL 4 treatments (not just the one it was selected from),
# so non-significant/absent signal in other treatments is visible too, not hidden.
# Two side-by-side panels (full analysis / cc-independent analysis), sharing the
# same pathway row order and the SAME color scale (unlike the DGE heatmap, FGSEA's
# NES is already normalized and comparable between the two frameworks).
#
# Usage: Rscript 11_cross_treatment_pathway_bubble.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
  library(dplyr)
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
source(file.path(proj_dir, "scripts", "utils_de_pathway.R"))

TOP_N_PER_TREATMENT <- 8

# `clean_pathway_name()` and `selection_rules` (the same per-treatment relevance
# rules as 09_fgsea_barplot_selected.R) are defined once in
# scripts/utils_de_pathway.R and shared here so the two plots' "curated"
# pathway sets cannot silently desync when one script is edited and not the
# other.

# --- Load every treatment's pathway_concordance files once (all pathways, any category) ---
load_treatment <- function(trt) {
  bind_rows(lapply(COLLECTIONS, function(coll) {
    fn <- file.path(comp_dir, paste0("pathway_concordance_", coll, "_Control_vs_", trt, ".csv"))
    read.csv(fn) %>% mutate(collection = coll)
  }))
}
all_data <- setNames(lapply(TREATMENTS, load_treatment), TREATMENTS)

# --- Per treatment: apply its own relevance rule, keep only the 3 significant/
#     concordant categories (matching 08/09; excludes sig_both_discordant and
#     not_sig_either), rank, take top N ---
pooled_ids <- list()
for (trt in TREATMENTS) {
  rule <- selection_rules[[trt]]
  inc_pat <- paste(rule$include, collapse = "|")
  exc_pat <- paste(rule$exclude, collapse = "|")
  cand <- all_data[[trt]] %>%
    filter(grepl(inc_pat, pathway, ignore.case = TRUE)) %>%
    filter(!grepl(exc_pat, pathway, ignore.case = TRUE)) %>%
    # Same 3-category allowlist as 08/09 (excludes sig_both_discordant), so this
    # plot's pooled pathway list can't silently include discordant-direction
    # pathways that the per-treatment barplots (08/09) would never show.
    filter(category %in% CATEGORIES) %>%
    mutate(best_padj = pmin(padj_full, padj_cc_independent, na.rm = TRUE)) %>%
    slice_min(best_padj, n = TOP_N_PER_TREATMENT, with_ties = FALSE)
  pooled_ids[[trt]] <- cand %>% distinct(collection, pathway)
  cat(sprintf("[%s] %d pathways selected for pooling\n", trt, nrow(pooled_ids[[trt]])))
}
pooled <- bind_rows(pooled_ids) %>% distinct(collection, pathway)
cat("\nPooled pathway list (union across treatments):", nrow(pooled), "pathways\n")

# --- Fetch this pooled list's values across ALL 4 treatments (incl. non-significant) ---
plot_data <- bind_rows(lapply(TREATMENTS, function(trt) {
  inner_join(pooled, all_data[[trt]], by = c("collection", "pathway")) %>%
    mutate(treatment = trt)
}))
plot_data$pathway_clean <- mapply(clean_pathway_name, plot_data$pathway, plot_data$collection)

# --- Row order: cluster by collection, then by mean full-analysis NES ---
row_order <- plot_data %>%
  group_by(collection, pathway_clean) %>%
  summarise(mean_nes = mean(NES_full, na.rm = TRUE), .groups = "drop") %>%
  arrange(collection, mean_nes) %>%
  pull(pathway_clean) %>%
  unique()

plot_data <- plot_data %>%
  mutate(
    pathway_clean = factor(pathway_clean, levels = row_order),
    treatment = factor(treatment, levels = TREATMENTS)
  )

nes_range <- range(c(plot_data$NES_full, plot_data$NES_cc_independent), na.rm = TRUE)
nes_limit <- max(abs(nes_range))

base_theme <- theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

make_panel <- function(nes_col, padj_col, title) {
  d <- plot_data
  d$NES <- d[[nes_col]]
  d$padj <- d[[padj_col]]
  d$neglog10padj <- pmin(-log10(pmax(d$padj, 1e-20)), 15)  # cap for display
  ggplot(d, aes(x = treatment, y = pathway_clean, size = neglog10padj, color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0,
                           limits = c(-nes_limit, nes_limit), name = "NES") +
    scale_size_continuous(range = c(0.5, 6), limits = c(0, 15), name = "-log10(padj)\n(capped at 15)") +
    labs(title = title, x = NULL, y = NULL) +
    base_theme
}

p_full <- make_panel("NES_full", "padj_full", "Full analysis")
p_cc   <- make_panel("NES_cc_independent", "padj_cc_independent", "Cell-cycle-independent analysis")

combined <- (p_full | p_cc) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Curated pathways across all 4 treatments",
    subtitle = paste(strwrap(paste0(
      "Top ", TOP_N_PER_TREATMENT, " biologically-relevant pathways per treatment (same relevance rules as fgsea_barplot_*_selected.png), pooled (", nrow(pooled), " unique pathways) and shown across ALL 4 treatments — ",
      "small pale dots mean that pathway simply isn't significant in that treatment, not that it wasn't tested. NES color scale is shared between panels (comparable across full/cc-independent)."
    ), width = 160), collapse = "\n"),
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) & theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "bottom")

out_file <- file.path(out_dir, "cross_treatment_pathway_bubble.png")
ggsave(out_file, combined, width = 14, height = max(8, nrow(pooled) * 0.28), dpi = 150, limitsize = FALSE, bg = "white")
cat("Saved:", out_file, "\n")

# ---------------------------------------------------------------------------
# Curated companion (2026-09-14, redesigned 2026-09-14) — same significance
# gate as above, but pathways selected by biological theme (phenotype_themes
# / select_curated_pathways() in scripts/utils_de_pathway.R, applied
# identically across treatments) and ranked by |NES|, instead of the
# per-treatment keyword `selection_rules` + best_padj ranking used above.
# Single-panel layout (shared left theme-label + pathway-name columns, no
# faceting) via save_curated_pathway_bubble_figure() in
# scripts/utils_de_pathway.R. Additional plot only ("_curated" suffix) —
# everything above (including `pooled`/`plot_data`/`out_file`) is untouched.
# ---------------------------------------------------------------------------
pooled_ids_curated <- list()
for (trt in TREATMENTS) {
  cand <- select_curated_pathways(all_data[[trt]])
  pooled_ids_curated[[trt]] <- cand %>% distinct(collection, pathway, theme, themes)
  cat(sprintf("[%s] %d themed pathways selected for pooling (curated)\n", trt, nrow(pooled_ids_curated[[trt]])))
}
# A pathway's theme membership is name-based and so agrees across treatments,
# but which theme it's FILED under (when it matches >1) is picked from
# whichever treatment's selection is encountered first, in TREATMENTS order
# (IL, TGF, PDGF_B, PDGF_C) — deterministic and reproducible, though in the
# rare case a pathway's top-N-per-theme rank differs by treatment, its listed
# `theme` here reflects that first treatment's pick, not every treatment's.
pooled_curated <- bind_rows(pooled_ids_curated) %>%
  distinct(collection, pathway, .keep_all = TRUE)
cat("\nPooled curated pathway list (union across treatments):", nrow(pooled_curated), "pathways\n")

plot_data_curated <- bind_rows(lapply(TREATMENTS, function(trt) {
  inner_join(pooled_curated %>% select(collection, pathway, theme, themes),
             all_data[[trt]], by = c("collection", "pathway")) %>%
    mutate(treatment = trt)
}))

out_file_curated <- file.path(out_dir, "cross_treatment_pathway_bubble_curated.png")
save_curated_pathway_bubble_figure(
  plot_data_curated, TREATMENTS,
  title = "Curated (CAF-phenotype) pathways across all 4 treatments",
  out_file = out_file_curated
)
cat("Saved:", out_file_curated, "\n")
