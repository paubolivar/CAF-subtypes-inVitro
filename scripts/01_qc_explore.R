#!/usr/bin/env Rscript
# QC exploration only — no filtering applied here.
# Loads all 15 CAF2 Flex per-sample filtered matrices, computes standard
# per-cell QC metrics, and writes a summary table + diagnostic plots.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
data_dir <- file.path(proj_dir, "data", "CAF2_results", "per_sample_outs")
out_dir  <- file.path(proj_dir, "results", "qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
samples <- sort(samples)
cat("Found", length(samples), "samples:\n")
print(samples)

parse_meta <- function(s) {
  replicate <- regmatches(s, regexpr("r[0-9]+", s))
  condition <- sub("^CAF2_r[0-9]+_", "", s)
  data.frame(sample = s, replicate = replicate, condition = condition)
}

obj_list <- list()
for (s in samples) {
  h5_path <- file.path(data_dir, s, "sample_filtered_feature_bc_matrix.h5")
  cat("Loading", s, "...\n")
  mat <- Read10X_h5(h5_path)
  so <- CreateSeuratObject(counts = mat, project = s, min.cells = 0, min.features = 0)
  meta <- parse_meta(s)
  so$sample    <- meta$sample
  so$replicate <- meta$replicate
  so$condition <- meta$condition
  so[["percent.mt"]]   <- PercentageFeatureSet(so, pattern = "^MT-")
  so[["percent.ribo"]] <- PercentageFeatureSet(so, pattern = "^RP[SL]")
  obj_list[[s]] <- so
}

n_mt_genes   <- sum(grepl("^MT-",  rownames(obj_list[[1]])))
n_ribo_genes <- sum(grepl("^RP[SL]", rownames(obj_list[[1]])))
n_total_genes <- nrow(obj_list[[1]])
cat("\nGene panel size (per sample, identical panel):", n_total_genes, "\n")
cat("Mitochondrial genes in panel (^MT-):", n_mt_genes, "\n")
cat("Ribosomal protein genes in panel (^RP[SL]):", n_ribo_genes, "\n")

merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = samples)
merged$condition <- factor(merged$condition, levels = c("Control","IL","TGF","PDGF_B","PDGF_C"))
merged$replicate <- factor(merged$replicate, levels = c("r1","r2","r3"))

cat("\nTotal cells across all 15 samples (no filtering applied):", ncol(merged), "\n")

saveRDS(merged, file.path(proj_dir, "results", "caf2_merged_unfiltered.rds"))

# --- Per-sample summary table ---
qc_df <- merged@meta.data %>%
  select(sample, condition, replicate, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo)

summary_tbl <- qc_df %>%
  group_by(sample, condition, replicate) %>%
  summarise(
    n_cells        = n(),
    median_nCount  = median(nCount_RNA),
    mean_nCount    = round(mean(nCount_RNA), 1),
    sd_nCount      = round(sd(nCount_RNA), 1),
    p99_nCount     = round(quantile(nCount_RNA, 0.99), 1),
    median_nFeature= median(nFeature_RNA),
    mean_nFeature  = round(mean(nFeature_RNA), 1),
    median_pct_mt  = round(median(percent.mt), 3),
    mean_pct_mt    = round(mean(percent.mt), 3),
    max_pct_mt     = round(max(percent.mt), 3),
    median_pct_ribo= round(median(percent.ribo), 2),
    mean_pct_ribo  = round(mean(percent.ribo), 2),
    .groups = "drop"
  ) %>%
  arrange(condition, replicate)

write.csv(summary_tbl, file.path(out_dir, "qc_summary_per_sample.csv"), row.names = FALSE)
cat("\n--- Per-sample QC summary ---\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# --- Outlier flagging (descriptive only — MAD-based, per sample; NOT applied as filter) ---
mad_flag <- qc_df %>%
  group_by(sample) %>%
  mutate(
    nCount_med = median(nCount_RNA), nCount_mad = mad(nCount_RNA),
    nFeature_med = median(nFeature_RNA), nFeature_mad = mad(nFeature_RNA),
    nCount_high_outlier = nCount_RNA > (nCount_med + 5 * nCount_mad),
    nCount_low_outlier  = nCount_RNA < (nCount_med - 5 * nCount_mad),
    high_pct_mt = percent.mt > 5
  ) %>%
  ungroup()

outlier_summary <- mad_flag %>%
  group_by(sample, condition, replicate) %>%
  summarise(
    n_cells = n(),
    n_nCount_high_outlier_5mad = sum(nCount_high_outlier),
    pct_nCount_high_outlier_5mad = round(100 * mean(nCount_high_outlier), 2),
    n_nCount_low_outlier_5mad = sum(nCount_low_outlier),
    n_pct_mt_over5 = sum(high_pct_mt),
    pct_pct_mt_over5 = round(100 * mean(high_pct_mt), 3),
    .groups = "drop"
  ) %>%
  arrange(condition, replicate)

write.csv(outlier_summary, file.path(out_dir, "qc_outlier_flagging_descriptive.csv"), row.names = FALSE)
cat("\n--- Descriptive outlier flagging (5-MAD high nCount, mt%>5) — NOT filtered, informational only ---\n")
print(as.data.frame(outlier_summary), row.names = FALSE)

# --- Plots ---
p1 <- ggplot(qc_df, aes(x = sample, y = nCount_RNA, fill = condition)) +
  geom_violin(scale = "width") + scale_y_log10() +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "nCount_RNA per sample (log10 scale)", x = NULL, y = "UMI count")

p2 <- ggplot(qc_df, aes(x = sample, y = nFeature_RNA, fill = condition)) +
  geom_violin(scale = "width") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "nFeature_RNA per sample", x = NULL, y = "Genes detected")

p3 <- ggplot(qc_df, aes(x = sample, y = percent.mt, fill = condition)) +
  geom_violin(scale = "width") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "percent.mt per sample", x = NULL, y = "% mitochondrial")

p4 <- ggplot(qc_df, aes(x = sample, y = percent.ribo, fill = condition)) +
  geom_violin(scale = "width") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "percent.ribo per sample", x = NULL, y = "% ribosomal protein genes")

ggsave(file.path(out_dir, "qc_violins.png"), (p1 / p2 / p3 / p4), width = 12, height = 16, dpi = 150)

p5 <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA, color = condition)) +
  geom_point(size = 0.3, alpha = 0.3) + scale_x_log10() +
  facet_wrap(~sample, ncol = 5) + theme_minimal() +
  labs(title = "nFeature vs nCount per sample")
ggsave(file.path(out_dir, "qc_scatter_count_vs_feature.png"), p5, width = 16, height = 10, dpi = 150)

cells_per_sample <- qc_df %>% count(sample, condition, replicate)
p6 <- ggplot(cells_per_sample, aes(x = sample, y = n, fill = condition)) +
  geom_col() + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
  labs(title = "Cells passing Cell Ranger's own filtered-matrix call, per sample", y = "n cells")
ggsave(file.path(out_dir, "qc_cells_per_sample.png"), p6, width = 8, height = 5, dpi = 150)

cat("\nDone. Outputs written to:", out_dir, "\n")
cat("Merged unfiltered Seurat object saved to: results/caf2_merged_unfiltered.rds\n")
