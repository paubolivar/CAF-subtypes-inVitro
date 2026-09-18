#!/usr/bin/env Rscript
# Cross-signature summary (mean z-score version) — same as 17_cross_signature_bubble.R
# but built from the mean-z-score scoring (16_score_published_caf_signatures_zscore.R)
# instead of AddModuleScore, for direct side-by-side comparison of the two methods.
#
# Usage: Rscript 18_cross_signature_bubble_zscore.R [CELL_LINE]   (default: CAF2)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
args <- commandArgs(trailingOnly = TRUE)
CELL_LINE <- if (length(args) >= 1) args[1] else "CAF2"
cat("Cell line:", CELL_LINE, "\n")

sig_dir  <- file.path(proj_dir, "results", CELL_LINE, "caf_signature_scoring_zscore")

TREATMENTS <- c("IL", "TGF", "PDGF_B", "PDGF_C")
SIG_ALPHA <- 0.05

summary_df <- read.csv(file.path(sig_dir, "signature_treatment_summary.csv"))

# --- Row order: hierarchical clustering on each signature's 4-treatment delta
# vector (not grouped by source anymore) — groups signatures by how they actually
# respond, regardless of which paper they came from. Computed here (on z-score
# deltas, the more cross-signature-comparable metric) and saved so
# 17_cross_signature_bubble.R (AddModuleScore version) can reuse the identical
# order, keeping the two bubble plots directly comparable row-for-row.
#
# NOTE: this row order is exported purely for shared VISUALIZATION layout.
# 17_cross_signature_bubble.R's AddModuleScore statistics are an independent
# analysis and do not depend on this script's z-score results in any
# statistical sense -- only on which row each signature is plotted in.
delta_mat <- summary_df %>%
  select(signature, treatment, delta) %>%
  tidyr::pivot_wider(names_from = treatment, values_from = delta) %>%
  tibble::column_to_rownames("signature") %>%
  as.matrix()
hc <- hclust(dist(delta_mat), method = "average")
row_order <- rownames(delta_mat)[hc$order]
write.csv(data.frame(signature = row_order, order = seq_along(row_order)),
          file.path(sig_dir, "signature_cluster_order.csv"), row.names = FALSE)

plot_df <- summary_df %>%
  mutate(
    signature = factor(signature, levels = row_order),
    treatment = factor(treatment, levels = TREATMENTS),
    neglog10padj = pmin(-log10(pmax(padj, 1e-300)), 30),
    sig_flag = padj < SIG_ALPHA
  )

delta_limit <- max(abs(plot_df$delta))

p <- ggplot(plot_df, aes(x = treatment, y = signature, size = neglog10padj, color = delta)) +
  geom_point() +
  scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0,
                         limits = c(-delta_limit, delta_limit), name = "Mean z-score\ndelta vs Control") +
  scale_size_continuous(range = c(0.5, 6), name = "-log10(padj)\n(BH, capped at 30)") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 6.5),
    axis.text.x = element_text(size = 9, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Published CAF signature scores (mean z-score) vs Control, all treatments",
    subtitle = paste(strwrap(paste0(
      "72 signatures (Cords/Elyada/Gao/Mechta/Pietras/Wu/KPlab_curated; Ruth/Ostman/Ohlund excluded), ",
      "single-cell mean per-gene z-score aggregated to per-replicate means (15 replicate-level observations), tested with limma (moderated t, n=3/group vs Control) — replicate is the unit of inference, matching this project's pseudobulk DGE convention. ",
      "Rows hierarchically clustered by response pattern across treatments (not by source). Compare against cross_signature_bubble.png (AddModuleScore version, same row order)."
    ), width = 150), collapse = "\n"),
    x = NULL, y = NULL
  )

out_file <- file.path(sig_dir, "cross_signature_bubble_zscore.png")
ggsave(out_file, p, width = 11, height = max(10, length(row_order) * 0.16), dpi = 150, limitsize = FALSE, bg = "white")
cat("Saved:", out_file, "\n")

# --- Top 15 signatures by strongest single-treatment effect, for a quick-glance table ---
top_hits <- plot_df %>%
  filter(sig_flag) %>%
  arrange(desc(abs(delta))) %>%
  select(signature, source, treatment, mean_control, mean_treatment, delta, padj) %>%
  head(15)
cat("\n--- Top 15 signature x treatment effects (by |delta|, padj<0.05) ---\n")
print(top_hits, row.names = FALSE)
write.csv(top_hits, file.path(sig_dir, "cross_signature_top_hits.csv"), row.names = FALSE)
