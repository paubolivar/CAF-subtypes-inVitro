#!/usr/bin/env Rscript
# Combine published CAF gene signature files (different labs, different formats)
# into one tidy long-format CSV: source, signature, gene, rank.
# (logFC is intentionally NOT included in the output — only Wu's source file has
# it, and the user decided it's not needed/consistent enough to expose as a column.
# It's still used internally to compute Wu's rank, just not carried into the CSV.)
#
# Sources combined here:
#   - CAF_signatures_all_labs_16June_2020.txt (template format; Mechta/Ohlund/
#     Pietras/Ostman/Ruth signatures) — wide matrix, one column per signature,
#     genes listed down each column in rank order (implicit rank, no logFC).
#     File uses old-Mac \r-only line endings.
#   - Cords_CAFs.txt — same wide-matrix format, \r\n line endings, no logFC.
#   - Gao_CAFs.txt — same wide-matrix format, \r\n line endings, no logFC;
#     column names lack a lab prefix, so "Gao_" is prepended here.
#   - elyada_CAFs.csv — long format (cluster;gene, semicolon-delimited despite the
#     .csv extension, UTF-8 BOM, \r\n line endings), no logFC. The source file has
#     15 clusters spanning the whole tumor atlas (immune/epithelial/endothelial
#     included); per user decision, only the 6 CAF/stroma-relevant clusters are
#     kept here: iCAF, apCAF, myCAF, fibroblasts, perivascular, EMT-like.
#   - Wu_CAFs.txt — a Seurat FindAllMarkers()-style table (p_val, avg_logFC, pct.1,
#     pct.2, p_val_adj, cluster, gene; \r\n line endings), 4 clusters: myCAFs,
#     iCAFs, dPVL cells, imPVL cells (all CAF/perivascular-relevant, none dropped).
#     This IS the source with explicit logFC — rank is computed here by sorting
#     each cluster's genes by descending avg_logFC (highest logFC = rank 1), NOT
#     taken from file order (the file is already sorted per-cluster, but we
#     recompute explicitly rather than assume).
#
# For signatures without logFC, `rank` reflects the gene's position in the source
# file (its original, presumably-curated order).

suppressPackageStartupMessages({
  library(dplyr)
})

proj_dir <- "/Users/paulina/Documents/projects_bioinformatics/CAF-subtypes-inVitro"
sig_dir  <- file.path(proj_dir, "data", "published_data", "CAF_gene_signatures")

# --- Helper: read a "wide" gene-signature matrix regardless of line-ending style ---
read_wide_signature_matrix <- function(path) {
  raw <- readChar(path, file.info(path)$size, useBytes = TRUE)
  raw <- gsub("\r\n", "\n", raw)
  raw <- gsub("\r", "\n", raw)
  con <- textConnection(raw)
  df <- read.delim(con, sep = "\t", header = TRUE, fill = TRUE, quote = "",
                    stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
  close(con)
  df
}

# --- Helper: melt a wide gene-signature matrix into long format ---
melt_wide_signatures <- function(df, prefix = NULL) {
  bind_rows(lapply(names(df), function(col) {
    genes <- df[[col]]
    genes <- genes[!is.na(genes) & trimws(genes) != ""]
    if (length(genes) == 0) return(NULL)
    sig_name <- if (!is.null(prefix)) paste0(prefix, col) else col
    data.frame(
      signature = sig_name,
      gene = trimws(genes),
      rank = seq_along(genes),
      logFC = NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

# ---------------------------------------------------------------------------
# 1. Template: Mechta / Ohlund / Pietras / Ostman / Ruth
# ---------------------------------------------------------------------------
template_df <- read_wide_signature_matrix(file.path(sig_dir, "CAF_signatures_all_labs_16June_2020.txt"))
template_long <- melt_wide_signatures(template_df)
cat("Template:", ncol(template_df), "signatures,", nrow(template_long), "gene rows\n")

# ---------------------------------------------------------------------------
# 2. Cords
# ---------------------------------------------------------------------------
cords_df <- read_wide_signature_matrix(file.path(sig_dir, "Cords_CAFs.txt"))
cords_long <- melt_wide_signatures(cords_df)
cat("Cords:", ncol(cords_df), "signatures,", nrow(cords_long), "gene rows\n")

# ---------------------------------------------------------------------------
# 3. Gao (column names lack a lab prefix -> prepend "Gao_")
# ---------------------------------------------------------------------------
gao_df <- read_wide_signature_matrix(file.path(sig_dir, "Gao_CAFs.txt"))
gao_long <- melt_wide_signatures(gao_df, prefix = "Gao_")
cat("Gao:", ncol(gao_df), "signatures,", nrow(gao_long), "gene rows\n")

# ---------------------------------------------------------------------------
# 4. Elyada (long format; keep only the 6 CAF/stroma-relevant clusters)
# ---------------------------------------------------------------------------
elyada_raw <- readChar(file.path(sig_dir, "elyada_CAFs.csv"),
                        file.info(file.path(sig_dir, "elyada_CAFs.csv"))$size, useBytes = TRUE)
elyada_raw <- sub("^\xef\xbb\xbf", "", elyada_raw, useBytes = TRUE)  # strip UTF-8 BOM
elyada_raw <- gsub("\r\n", "\n", elyada_raw)
elyada_con <- textConnection(elyada_raw)
elyada_df <- read.delim(elyada_con, sep = ";", header = TRUE, quote = "", stringsAsFactors = FALSE)
close(elyada_con)

ELYADA_KEEP <- c("iCAF", "apCAF", "myCAF", "fibroblasts", "perivascular", "EMT-like")
elyada_long <- elyada_df %>%
  filter(cluster %in% ELYADA_KEEP) %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  transmute(signature = paste0("Elyada_", cluster), gene = trimws(Associated.Gene.Name), rank, logFC = NA_real_)
cat("Elyada:", length(ELYADA_KEEP), "signatures kept (of", length(unique(elyada_df$cluster)), "total in source),",
    nrow(elyada_long), "gene rows\n")

# ---------------------------------------------------------------------------
# 5. Wu (has explicit logFC -> rank by descending avg_logFC per cluster)
# ---------------------------------------------------------------------------
wu_path <- file.path(sig_dir, "Wu_CAFs.txt")
wu_raw <- readChar(wu_path, file.info(wu_path)$size, useBytes = TRUE)
wu_raw <- gsub("\r\n", "\n", wu_raw)
wu_con <- textConnection(wu_raw)
wu_df <- read.delim(wu_con, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
close(wu_con)

wu_long <- wu_df %>%
  group_by(cluster) %>%
  arrange(desc(avg_logFC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  transmute(
    signature = paste0("Wu_", gsub(" ", "_", cluster)),
    gene = trimws(gene),
    rank,
    logFC = avg_logFC
  )
cat("Wu:", length(unique(wu_df$cluster)), "signatures,", nrow(wu_long), "gene rows (ranked by descending avg_logFC)\n")

# ---------------------------------------------------------------------------
# 6. KPlab_curated — hand-curated by the user's colleague, no logFC/rank given;
#    rank = position as given (list order), already deduplicated (a few input
#    lists had repeated genes — see KPlab_curated_CAFs.csv comments/chat context).
#    Signature names already carry the "KPlab_curated_" prefix in the source file
#    since "KPlab_curated" itself contains an underscore, which would break the
#    generic source = sub("_.*","",signature) rule below — handled explicitly.
# ---------------------------------------------------------------------------
kplab_df <- read.csv(file.path(sig_dir, "KPlab_curated_CAFs.csv"), stringsAsFactors = FALSE)
kplab_long <- kplab_df %>%
  group_by(signature) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(logFC = NA_real_) %>%
  select(signature, gene, rank, logFC)
cat("KPlab_curated:", length(unique(kplab_long$signature)), "signatures,", nrow(kplab_long), "gene rows\n")

# ---------------------------------------------------------------------------
# Combine + derive `source` (lab) from the signature name prefix
# ---------------------------------------------------------------------------
TOP_N_GENES <- 30  # each signature is capped to its top N genes (by rank); user
                    # decision 2026-09-11 based on prior experience with CAF
                    # signatures — was 50, adjustable for future exploration

combined <- bind_rows(template_long, cords_long, gao_long, elyada_long, wu_long, kplab_long) %>%
  mutate(source = ifelse(grepl("^KPlab_curated_", signature), "KPlab_curated", sub("_.*", "", signature))) %>%
  filter(rank <= TOP_N_GENES) %>%
  select(source, signature, gene, rank) %>%
  arrange(source, signature, rank)

out_file <- file.path(sig_dir, "CAF_gene_signatures_combined.csv")
write.csv(combined, out_file, row.names = FALSE)

cat("\n--- Summary ---\n")
cat("Total rows:", nrow(combined), "\n")
cat("Total unique signatures:", length(unique(combined$signature)), "\n")
print(combined %>% distinct(source, signature) %>% count(source, name = "n_signatures"))
cat("\nSaved:", out_file, "\n")
