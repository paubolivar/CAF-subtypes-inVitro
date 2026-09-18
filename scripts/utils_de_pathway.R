# Shared helper functions for the DE / pathway-analysis pipeline (scripts
# 05, 06, 08, 09, 11). Source this after 00_pipeline_config.R. Consolidates
# logic that was previously copy-pasted near-verbatim across scripts, so a fix
# or a curated-pathway-list edit can no longer silently apply to one script
# but not the other.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(msigdbr)
  library(fgsea)
  library(patchwork)
})

# --- MSigDB collection resolution (05, 06) ----------------------------------
# KEGG's subcollection name has changed across msigdbr releases; resolve it
# dynamically instead of hardcoding a name that might not exist in the
# installed msigdbr version.
resolve_collections <- function() {
  coll_info <- msigdbr_collections()
  kegg_subcol <- unique(grep("^CP:KEGG", coll_info$gs_subcollection, value = TRUE))
  kegg_subcol <- if ("CP:KEGG_LEGACY" %in% kegg_subcol) "CP:KEGG_LEGACY" else kegg_subcol[1]
  cat("\nUsing KEGG subcollection:", kegg_subcol, "\n")
  list(
    Hallmark  = list(collection = "H",  subcollection = NULL),
    GO_BP     = list(collection = "C5", subcollection = "GO:BP"),
    Oncogenic = list(collection = "C6", subcollection = NULL),
    Reactome  = list(collection = "C2", subcollection = "CP:REACTOME"),
    KEGG      = list(collection = "C2", subcollection = kegg_subcol)
  )
}

get_pathways <- function(coll, subcoll) {
  gs <- msigdbr(species = "Homo sapiens", collection = coll, subcollection = subcoll)
  split(gs$gene_symbol, gs$gs_name)
}

# --- FGSEA runner (05, 06) ---------------------------------------------------
# Runs fgsea for every MSigDB collection x treatment, ranking genes by
# `stat_col` (DESeq2's Wald `stat` in 05, limma's moderated `t` in 06 -- both
# are "test statistics", the matched ranking metric for comparability between
# the two DE engines). Writes one CSV per collection/treatment plus a combined
# `fgsea_all.csv` to `out_dir`, and returns the combined data frame.
run_fgsea_all <- function(all_res, treatments, stat_col, out_dir) {
  collections <- resolve_collections()
  fgsea_all <- list()
  for (coll_name in names(collections)) {
    info <- collections[[coll_name]]
    pathways <- get_pathways(info$collection, info$subcollection)
    cat(sprintf("\n%s: %d gene sets\n", coll_name, length(pathways)))

    for (trt in treatments) {
      df <- all_res[[trt]]
      df <- df[!is.na(df[[stat_col]]), ]
      ranks <- sort(setNames(df[[stat_col]], df$gene), decreasing = TRUE)
      fres <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500, eps = 0)
      fres$leadingEdge <- vapply(fres$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
      fres$collection <- coll_name
      fres$contrast <- paste0("Control_vs_", trt)
      fres <- fres[order(fres$pval), ]
      write.csv(fres, file.path(out_dir, paste0("fgsea_", coll_name, "_Control_vs_", trt, ".csv")), row.names = FALSE)
      fgsea_all[[paste(coll_name, trt)]] <- fres
      n_sig <- sum(fres$padj < SIG_ALPHA, na.rm = TRUE)
      cat(sprintf("  %s vs Control: %d pathways significant at padj<%s\n", trt, n_sig, SIG_ALPHA))
    }
  }
  fgsea_all_df <- bind_rows(fgsea_all)
  write.csv(fgsea_all_df, file.path(out_dir, "fgsea_all.csv"), row.names = FALSE)
  fgsea_all_df
}

# --- FGSEA dotplots (05, 06) -------------------------------------------------
write_fgsea_dotplots <- function(fgsea_all_df, out_dir, title_suffix) {
  for (coll_name in unique(fgsea_all_df$collection)) {
    sub <- fgsea_all_df[fgsea_all_df$collection == coll_name, ]
    top <- sub %>%
      filter(!is.na(padj), padj < SIG_ALPHA) %>%
      group_by(contrast) %>%
      slice_max(order_by = abs(NES), n = PATHWAY_TOP_N, with_ties = FALSE) %>%
      ungroup()
    if (nrow(top) == 0) next
    top$pathway <- factor(top$pathway, levels = rev(unique(top$pathway)))
    p <- ggplot(top, aes(x = NES, y = pathway, fill = NES > 0)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      facet_wrap(~contrast, scales = "free_y", ncol = 2) +
      theme_minimal(base_size = 8) +
      labs(title = paste0("Top FGSEA pathways — ", coll_name, " (", title_suffix, ")"),
           x = "Normalized enrichment score (NES)", y = NULL) +
      theme(legend.position = "none")
    ggsave(file.path(out_dir, paste0("fgsea_dotplot_", coll_name, ".png")), p,
           width = 14, height = 10, dpi = 150, limitsize = FALSE)
  }
}

# --- Volcano plot (05, 06) ---------------------------------------------------
# `df` must already carry a `sig` logical column -- the significance/effect-
# size cutoff differs by DE engine (see fc_cutoff() in 00_pipeline_config.R)
# and stays in the calling script; this shares only the plotting and
# top-12-label-selection logic, which is identical between engines: label the
# top 12 significant genes, sorted by ascending `padj_col`, tie-broken by
# descending |effect size|.
build_volcano_plot <- function(df, x_col, y_pval_col, padj_col, title, xlab) {
  df$.x <- df[[x_col]]
  df$.y <- -log10(df[[y_pval_col]])
  df$.padj <- df[[padj_col]]

  label_df <- df %>%
    filter(sig) %>%
    arrange(.padj, desc(abs(.x))) %>%
    slice_head(n = 12)

  ggplot(df, aes(x = .x, y = .y, color = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_text_repel(data = label_df, aes(label = gene),
                    size = 2.5, max.overlaps = 20, box.padding = 0.35,
                    point.padding = 0.2, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    theme_minimal() +
    labs(title = title, x = xlab, y = "-log10(p)") +
    theme(legend.position = "none")
}

# --- Pathway concordance category display (08, 09, 11) ----------------------
CATEGORIES <- c("sig_full_only", "sig_cc_independent_only", "sig_both_concordant")

cat_colors <- c(
  sig_full_only           = "steelblue",
  sig_cc_independent_only = "lightcoral",
  sig_both_concordant     = "mediumpurple"
)
cat_labels <- c(
  sig_full_only           = "Full analysis only (CC included)",
  sig_cc_independent_only = "CC-independent analysis only",
  sig_both_concordant     = "Significant in both (robust)"
)

clean_pathway_name <- function(pw, coll) {
  pw <- switch(coll,
    Hallmark = sub("^HALLMARK_", "", pw),
    GO_BP    = sub("^GOBP_", "", pw),
    Reactome = sub("^REACTOME_", "", pw),
    KEGG     = sub("^KEGG_", "", pw),
    pw
  )
  gsub("_", " ", pw)
}

# --- Per-treatment biological relevance rule for curated pathway selection
# (09_fgsea_barplot_selected.R and 11_cross_treatment_pathway_bubble.R) --
# single source of truth so the "curated" pathway set shown in the
# per-treatment barplots and the cross-treatment bubble plot cannot silently
# desync from one another when one script is edited and not the other.
selection_rules <- list(
  IL = list(
    include = c("INTERFERON", "INFLAMMAT", "NFKB", "NF_KB", "TNFA", "^TNF_", "_TNF_",
                "CYTOKINE", "INTERLEUKIN", "\\bIL[0-9]", "CHEMOKINE", "COMPLEMENT",
                "\\bJAK\\b", "\\bSTAT[0-9]?\\b", "ALLOGRAFT", "ANTIGEN", "T_HELPER",
                "TH17", "T_CELL", "LYMPHOCYTE", "\\bIMMUNE\\b", "TOLL_LIKE",
                "NOD_LIKE", "RIG_I", "MYD88"),
    exclude = c("VIRUS", "VIRAL", "PRION", "BACTERI"),
    include_readable = "interferon, inflammation, NFkB, TNF, cytokine, interleukin, chemokine, complement, JAK/STAT, allograft rejection, antigen presentation, T-cell/T-helper/lymphocyte, general immune, Toll-like/NOD-like/RIG-I receptor signaling",
    exclude_readable = "generic viral-defense-labeled, prion-disease, and bacterial-infection gene sets"
  ),
  TGF = list(
    include = c("TGF_BETA", "\\bTGFB", "EPITHELIAL_MESENCHYMAL", "\\bEMT\\b",
                "COLLAGEN", "EXTRACELLULAR_MATRIX", "\\bECM\\b", "\\bMATRIX\\b",
                "PROTEOGLYCAN", "MYOFIBROBLAST", "FIBROSIS", "FIBROBLAST",
                "SYNDECAN", "UNFOLDED_PROTEIN", "\\bUPR\\b", "INTEGRIN",
                "SMOOTH_MUSCLE", "MYOGENESIS", "FOCAL_ADHESION", "WOUND",
                "\\bCSR_", "INTERFERON"),
    exclude = c("VIRUS", "VIRAL", "PRION", "BACTERI"),
    include_readable = "TGF-beta targets, EMT, collagen/ECM/proteoglycans, myofibroblast/fibrosis/fibroblast, syndecan, unfolded protein response, integrin/focal adhesion, smooth muscle/myogenesis, wound healing, Core Serum Response (fibroblast activation signature), interferon signaling (to show the suppression finding and its robustness)",
    exclude_readable = "generic viral-defense-labeled, prion-disease, and bacterial-infection gene sets"
  ),
  PDGF_B = list(
    include = c("COLLAGEN", "EXTRACELLULAR_MATRIX", "\\bECM\\b", "\\bMATRIX\\b",
                "PROTEOGLYCAN", "SMOOTH_MUSCLE", "MYOGENESIS", "FOCAL_ADHESION",
                "CELL_CYCLE", "MITOTIC", "\\bE2F\\b", "\\bMYC\\b", "DNA_REPLICATION",
                "RIBOSOME", "\\bRRNA\\b", "CHROMATID", "\\bKRAS\\b", "\\bRAS\\.",
                "MAPK", "PROTEASOME", "\\bCSR_", "EPITHELIAL_MESENCHYMAL"),
    exclude = c("VIRUS", "VIRAL", "PRION", "BACTERI"),
    include_readable = "ECM/collagen/proteoglycans, smooth muscle/myogenesis, focal adhesion (the robust suppression side); cell cycle/mitosis/E2F/MYC/DNA replication/ribosome biogenesis (the proliferation side, expect mostly cell-cycle-confounded); KRAS/RAS/MAPK signaling, proteasome, Core Serum Response, EMT",
    exclude_readable = "generic viral-defense-labeled, prion-disease, and bacterial-infection gene sets"
  ),
  PDGF_C = list(
    include = c("COLLAGEN", "EXTRACELLULAR_MATRIX", "\\bECM\\b", "\\bMATRIX\\b",
                "PROTEOGLYCAN", "SMOOTH_MUSCLE", "MYOGENESIS", "FOCAL_ADHESION",
                "CELL_CYCLE", "MITOTIC", "\\bE2F\\b", "\\bMYC\\b", "DNA_REPLICATION",
                "RIBOSOME", "\\bRRNA\\b", "CHROMATID", "\\bKRAS\\b", "\\bRAS\\.",
                "MAPK", "PROTEASOME", "\\bCSR_", "EPITHELIAL_MESENCHYMAL"),
    exclude = c("VIRUS", "VIRAL", "PRION", "BACTERI"),
    include_readable = "same rule as PDGF_B (ECM/collagen/matrix, smooth muscle/myogenesis, focal adhesion, cell cycle/proliferation machinery, KRAS/RAS/MAPK, proteasome, Core Serum Response, EMT) — kept identical so PDGF_B vs PDGF_C potency is directly comparable",
    exclude_readable = "generic viral-defense-labeled, prion-disease, and bacterial-infection gene sets"
  )
)

# --- Biologically curated pathway themes ("_curated" FGSEA plots, 2026-09-14) ---
# Treatment-AGNOSTIC (unlike `selection_rules` above, which varies per treatment):
# applied identically to every treatment/comparison, tied to the project's core
# biological question — do IL/TNF/PDGF-family treatments modify CAF phenotype/
# state, and in what direction? Regex uses biological synonyms/annotations, not
# literal treatment names. Catalog order is the deterministic tie-break for a
# pathway's *primary* theme when its name matches more than one (see
# select_curated_pathways() below).
phenotype_themes <- list(
  Interferon_Response = list(
    include = c("INTERFERON", "\\bIFN\\b", "\\bISG\\b", "\\bIRF[0-9]?\\b",
                "\\bJAK\\b", "\\bSTAT[0-9]?\\b"),
    exclude = c("PRION", "VIRUS", "VIRAL", "BACTERI")
  ),
  Inflammatory_Response = list(
    include = c("INFLAMMAT", "ACUTE_PHASE", "IMMUNE_RESPONSE", "LEUKOCYTE",
                "COMPLEMENT", "\\bIMMUNE\\b"),
    exclude = c("PRION", "VIRUS", "VIRAL", "BACTERI")
  ),
  Cytokine_Signaling = list(
    include = c("CYTOKINE", "INTERLEUKIN", "\\bIL[0-9]", "CHEMOKINE"),
    exclude = character(0)
  ),
  TNF_NFkB = list(
    include = c("\\bTNF\\b", "TNFA", "NFKB", "NF_KB", "\\bRELA\\b", "\\bRELB\\b",
                "\\bIKK", "\\bNIK\\b"),
    exclude = character(0)
  ),
  TGFbeta_Signaling = list(
    include = c("TGF_BETA", "\\bTGFB", "\\bSMAD[0-9]?\\b", "BMP_SIGNALING"),
    exclude = character(0)
  ),
  ECM_Organization = list(
    include = c("EXTRACELLULAR_MATRIX", "\\bECM\\b", "MATRISOME",
                "BASEMENT_MEMBRANE", "EXTRACELLULAR_STRUCTURE"),
    exclude = character(0)
  ),
  Collagen = list(
    include = c("COLLAGEN", "FIBRILLOGENESIS", "PROTEOGLYCAN"),
    exclude = character(0)
  ),
  Fibroblast_Activation = list(
    include = c("FIBROBLAST", "FIBROSIS", "MYOFIBROBLAST", "STELLATE_CELL"),
    exclude = character(0)
  ),
  Cell_Adhesion = list(
    include = c("CELL_ADHESION", "FOCAL_ADHESION", "\\bINTEGRIN", "CADHERIN",
                "ADHESION_MOLECULE"),
    exclude = character(0)
  ),
  Wound_Healing_Remodeling = list(
    include = c("WOUND", "TISSUE_REMODELING", "REGENERATION"),
    exclude = character(0)
  ),
  Matrix_Degradation_Proteolysis = list(
    include = c("METALLOPROTEINASE", "\\bMMP[0-9]*\\b", "PROTEOLYSIS",
                "PEPTIDASE", "ECM_DISASSEMBLY", "COLLAGEN_CATABOLIC",
                "COLLAGEN_DEGRADATION"),
    exclude = character(0)
  ),
  PDGF_GrowthFactor_Signaling = list(
    include = c("\\bPDGF\\b", "PLATELET_DERIVED_GROWTH_FACTOR",
                "GROWTH_FACTOR_RECEPTOR", "GROWTH_FACTOR_BINDING",
                "RECEPTOR_TYROSINE_KINASE"),
    exclude = character(0)
  ),
  PI3K_AKT_MAPK_Signaling = list(
    include = c("\\bPI3K\\b", "\\bAKT\\b", "\\bMTOR\\b", "\\bMAPK\\b",
                "\\bERK[0-9]?\\b", "\\bRAS\\b", "\\bRAF\\b", "\\bKRAS\\b"),
    exclude = character(0)
  ),
  Angiogenesis_Vascular = list(
    include = c("ANGIOGENESIS", "VASCULAR", "\\bVEGF", "ENDOTHELI"),
    exclude = character(0)
  ),
  Myofibroblast_Contractile = list(
    include = c("SMOOTH_MUSCLE", "MYOGENESIS", "CONTRACT", "ACTIN_CYTOSKELETON",
                "\\bMYOSIN\\b"),
    exclude = character(0)
  ),
  EMT_Mesenchymal = list(
    include = c("EPITHELIAL_MESENCHYMAL", "\\bEMT\\b", "MESENCHYMAL"),
    exclude = character(0)
  ),
  Hypoxia_Metabolic = list(
    include = c("HYPOXIA", "GLYCOLYSIS", "OXIDATIVE_PHOSPHORYLATION",
                "\\bMETABOL"),
    exclude = character(0)
  ),
  # Added 2026-09-14, deliberately kept small (max_total below) -- captures
  # canonical proliferation/cell-cycle-machinery programs (E2F targets,
  # G2/M checkpoint, mitotic spindle/mitosis, DNA replication/S phase,
  # MYC-associated proliferation), NOT generic cell-type-specific
  # "proliferation" (that would collide with e.g. Angiogenesis_Vascular's
  # "endothelial cell proliferation" or Myofibroblast_Contractile's "smooth
  # muscle cell proliferation" hits -- deliberately no bare "PROLIFERATION"
  # keyword here for that reason), and explicitly not generic DNA repair,
  # chromatin remodeling, apoptosis, or bare p53 signaling.
  Proliferation_Cell_Cycle = list(
    include = c("\\bE2F\\b", "E2F[0-9]?_TARGETS", "G2M_CHECKPOINT", "\\bG2M\\b",
                "G2_M_", "MITOTIC_SPINDLE", "\\bMITOSIS\\b", "\\bMITOTIC\\b",
                "DNA_REPLICATION", "\\bS_PHASE\\b", "\\bCELL_CYCLE\\b",
                "MYC_TARGETS", "CHROMOSOME_SEGREGATION", "SPINDLE_ASSEMBLY",
                "KINETOCHORE", "SISTER_CHROMATID"),
    exclude = c("DNA_REPAIR", "MISMATCH_REPAIR", "NUCLEOTIDE_EXCISION_REPAIR",
                "BASE_EXCISION_REPAIR", "HOMOLOGOUS_RECOMBINATION_REPAIR",
                "APOPTOSIS", "APOPTOTIC", "P53_PATHWAY", "\\bTP53\\b",
                "CHROMATIN_REMODELING", "CHROMATIN_ORGANIZATION",
                "HISTONE_MODIFICATION"),
    max_total = 5, max_per_collection = 2
  )
)

# Selects, from an already-FDR-filtered pathway_concordance data frame (one
# treatment, all 5 collections stacked with a `collection` column — same shape
# as the `pathway_concordance_<coll>_Control_vs_<TRT>.csv` files read
# elsewhere), the pathways matching each `phenotype_themes` entry. FDR gate is
# unchanged (`category %in% CATEGORIES`, i.e. padj_full/padj_cc_independent <
# SIG_ALPHA, direction-concordant if significant in both — never raw p-value).
# Within a theme: capped at `top_n_per_theme_collection` per (theme,
# collection) first (redundancy control, so GO_BP's many near-duplicate terms
# can't crowd out the other 4 collections), then `top_n_per_theme` overall by
# |NES|. Either cap can be overridden per-theme via `max_total`/
# `max_per_collection` fields on that theme's entry in `phenotype_themes`
# (e.g. Proliferation_Cell_Cycle's max_total=5) — themes without those fields
# fall back to this function's `top_n_per_theme`/`top_n_per_theme_collection`
# defaults, unaffected. A pathway matching multiple themes' regex is kept
# once — placed in its first-matching theme in catalog order above
# (deterministic) — with every matched theme recorded in the `themes` column
# for transparency. Used by 08_fgsea_barplots_by_treatment.R and
# 11_cross_treatment_pathway_bubble.R to build their `_curated` companion
# plots.
select_curated_pathways <- function(df, themes = phenotype_themes,
                                     top_n_per_theme = 8,
                                     top_n_per_theme_collection = 3) {
  df <- df %>%
    filter(category %in% CATEGORIES) %>%
    mutate(display_NES = ifelse(category == "sig_cc_independent_only", NES_cc_independent, NES_full))

  hits_list <- list()
  for (theme_name in names(themes)) {
    rule <- themes[[theme_name]]
    inc_pat <- paste(rule$include, collapse = "|")
    hits <- df[grepl(inc_pat, df$pathway, ignore.case = TRUE), ]
    if (length(rule$exclude) > 0) {
      exc_pat <- paste(rule$exclude, collapse = "|")
      hits <- hits[!grepl(exc_pat, hits$pathway, ignore.case = TRUE), ]
    }
    if (nrow(hits) == 0) next
    this_top_n <- if (!is.null(rule$max_total)) rule$max_total else top_n_per_theme
    this_top_n_coll <- if (!is.null(rule$max_per_collection)) rule$max_per_collection else top_n_per_theme_collection
    hits <- hits %>%
      group_by(collection) %>%
      slice_max(abs(display_NES), n = this_top_n_coll, with_ties = FALSE) %>%
      ungroup() %>%
      slice_max(abs(display_NES), n = this_top_n, with_ties = FALSE)
    hits$theme <- theme_name
    hits_list[[theme_name]] <- hits
  }
  if (length(hits_list) == 0) return(df[0, ])
  all_hits <- bind_rows(hits_list)

  # A pathway's NAME-based theme membership is deterministic, but WHICH themes
  # it survives the top-N cap under can vary by input (i.e. by treatment, when
  # this is called once per treatment) since ranking is by that treatment's
  # |NES|. `themes` records every theme it matched in THIS call; `theme` (the
  # single facet it's plotted under) is the first of those in catalog order.
  membership <- all_hits %>%
    group_by(collection, pathway, contrast) %>%
    summarise(themes = paste(unique(theme), collapse = ", "), .groups = "drop")

  all_hits %>%
    distinct(collection, pathway, contrast, .keep_all = TRUE) %>%
    select(-theme) %>%
    left_join(membership, by = c("collection", "pathway", "contrast")) %>%
    mutate(theme = sub(",.*$", "", themes))
}

# --- Single-panel curated pathway figure (2026-09-14 redesign) --------------
# Replaces the earlier facet_wrap(~theme) layout (too many narrow panels) with
# ONE horizontal-bar panel plus a narrow left column of bold, vertically-
# centered biological-theme labels — visualization only. Consumes exactly the
# same select_curated_pathways() output; does not alter selection, FDR
# gating, ranking, or any NES/padj value. Used by
# 08_fgsea_barplots_by_treatment.R (per-treatment) and (adapted for a
# treatment x pathway dot grid) 11_cross_treatment_pathway_bubble.R.

# Fixed, biologically meaningful display order (NOT alphabetical) for the
# left-side theme column -- presentation only, independent of
# `phenotype_themes`' catalog order (which only controls the deterministic
# primary-theme tie-break in select_curated_pathways(), not display order). A
# theme absent from a given treatment's curated selection is simply skipped,
# never padded in as an empty section.
THEME_DISPLAY_ORDER <- c(
  "Interferon_Response", "Inflammatory_Response", "Cytokine_Signaling",
  "TNF_NFkB", "TGFbeta_Signaling", "Fibroblast_Activation",
  "Myofibroblast_Contractile", "ECM_Organization", "Collagen",
  "Matrix_Degradation_Proteolysis", "Cell_Adhesion",
  "Wound_Healing_Remodeling", "PDGF_GrowthFactor_Signaling",
  "PI3K_AKT_MAPK_Signaling", "Proliferation_Cell_Cycle", "EMT_Mesenchymal",
  "Angiogenesis_Vascular", "Hypoxia_Metabolic"
)
stopifnot(setequal(THEME_DISPLAY_ORDER, names(phenotype_themes)))

THEME_DISPLAY_LABELS <- c(
  Interferon_Response            = "Interferon\nresponse",
  Inflammatory_Response          = "Inflammatory\nresponse",
  Cytokine_Signaling             = "Cytokine\nsignaling",
  TNF_NFkB                       = "TNF /\nNF-kB",
  TGFbeta_Signaling              = "TGF-beta\nsignaling",
  Fibroblast_Activation          = "Fibroblast\nactivation",
  Myofibroblast_Contractile      = "Myofibroblast /\ncontractile",
  ECM_Organization               = "ECM\norganization",
  Collagen                       = "Collagen",
  Matrix_Degradation_Proteolysis = "Matrix degradation /\nproteolysis",
  Cell_Adhesion                  = "Cell\nadhesion",
  Wound_Healing_Remodeling       = "Wound healing /\ntissue remodeling",
  PDGF_GrowthFactor_Signaling    = "PDGF / growth-\nfactor signaling",
  PI3K_AKT_MAPK_Signaling        = "PI3K / AKT /\nMAPK",
  Proliferation_Cell_Cycle       = "Proliferation /\ncell cycle",
  EMT_Mesenchymal                = "EMT /\nmesenchymal",
  Angiogenesis_Vascular          = "Angiogenesis /\nvascular",
  Hypoxia_Metabolic              = "Hypoxia /\nmetabolic"
)

wrap_label <- function(x, width = 46) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

# Shared row-layout step: given a curated (or curated-derived) data frame with
# a `theme` column, returns it with `theme` releveled to the realized subset
# of THEME_DISPLAY_ORDER, plus row positions `y` (largest at the top) ordered
# by theme (fixed order) then by `order_col` (a numeric column name, e.g.
# "display_NES" or "mean_abs_nes") descending in |value|. A `gap` of extra
# vertical units is inserted BETWEEN consecutive theme groups -- genuine
# whitespace (not just a thin separator line), so a short (e.g.
# single-pathway) theme group doesn't visually crowd its neighbors, no matter
# how few rows it has.
.assign_theme_rows <- function(df, order_col, gap = 0.9) {
  realized_order <- THEME_DISPLAY_ORDER[THEME_DISPLAY_ORDER %in% unique(df$theme)]
  df <- df %>%
    mutate(theme = factor(theme, levels = realized_order)) %>%
    arrange(theme, desc(abs(.data[[order_col]]))) %>%
    mutate(row_id = row_number())
  theme_ord <- as.integer(df$theme)
  pos <- df$row_id + (theme_ord - 1) * gap
  df$y <- max(pos) - pos + min(pos)  # reverse so theme 1 / row 1 ends up at the top
  df
}

# Bold, vertically-centered theme label metadata + alternating background
# bands (the gap from .assign_theme_rows() already creates real whitespace
# between bands, so no separate boundary-line layer is needed). Shared by
# every panel in a figure -- all panels must use the identical `ylims` (and
# the identical upstream `y` column) to stay row-aligned.
.theme_group_geometry <- function(df) {
  theme_bounds <- df %>%
    group_by(theme) %>%
    summarise(y_min = min(y), y_max = max(y), y_mid = mean(y), .groups = "drop") %>%
    mutate(label = THEME_DISPLAY_LABELS[as.character(theme)],
           band = seq_len(n()) %% 2 == 0)
  pad <- 1.3
  list(theme_bounds = theme_bounds, ylims = c(min(df$y) - pad, max(df$y) + pad),
       band_rows = theme_bounds[theme_bounds$band, ])
}

.band_layer <- function(geo) {
  if (nrow(geo$band_rows) == 0) return(NULL)
  geom_rect(data = geo$band_rows, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = y_min - 0.5, ymax = y_max + 0.5),
            fill = "grey93")
}

.boundary_layer <- function(geo) {
  NULL
}

# A ggplot theme for a "label-only" left-column panel (theme names / pathway
# names) that reserves the SAME x-axis footprint as a real data panel --
# invisible (white) text, not element_blank(), so patchwork's row alignment
# between panels isn't thrown off by one side having no real x-axis/legend
# and the other having one. Pass the EXACT same axis_text_size/axis_title_size
# (and title_size, if the paired data panel has a plot title) used by the
# real data panel it's aligned with, so the reserved band heights match
# precisely even when fonts are large.
.blank_axis_theme <- function(base_size, axis_text_size = base_size,
                               axis_title_size = base_size + 1,
                               title_size = NULL) {
  base <- theme_minimal(base_size = base_size) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(size = axis_text_size, color = "white"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)
    )
  if (!is.null(title_size)) {
    base <- base + theme(axis.title.x = element_blank(),
                          plot.title = element_text(size = title_size, color = "white"))
  } else {
    base <- base + theme(axis.title.x = element_text(size = axis_title_size, color = "white"))
  }
  base
}

# --- 08_fgsea_barplots_by_treatment.R: one treatment's curated pathways ----
# Font sizes and dimensions tuned 2026-09-14 for readability at normal
# viewing size without zooming -- pure typography/layout, no change to
# selection, ranking, or values (see select_curated_pathways()).
save_curated_pathway_figure <- function(curated, title, out_file) {
  plot_df <- curated %>%
    mutate(pathway_label = wrap_label(paste0(mapply(clean_pathway_name, pathway, collection),
                                               "  [", collection, "]"), width = 50)) %>%
    .assign_theme_rows("display_NES")
  n <- nrow(plot_df)
  n_themes <- length(unique(plot_df$theme))
  geo <- .theme_group_geometry(plot_df)

  PATHWAY_LABEL_SIZE <- 12.5   # pt, axis.text.y
  THEME_LABEL_SIZE   <- 5.6    # mm, geom_text
  AXIS_TEXT_SIZE     <- 13     # pt
  AXIS_TITLE_SIZE    <- 15     # pt
  LEGEND_TEXT_SIZE   <- 13     # pt
  TITLE_SIZE         <- 18     # pt
  SUBTITLE_SIZE      <- 12.5   # pt

  p_bars <- ggplot(plot_df, aes(x = display_NES, y = y)) +
    .band_layer(geo) +
    # orientation="y" is required: with a numeric (not discrete-factor) `y`,
    # geom_col() cannot auto-detect which axis is the bar-grouping axis and
    # silently mis-stacks bars across rows (position_stack() warning) instead
    # of drawing one bar per row.
    geom_col(aes(fill = category), width = 0.55, orientation = "y", position = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    .boundary_layer(geo) +
    scale_y_continuous(breaks = plot_df$y, labels = plot_df$pathway_label, expand = c(0, 0)) +
    coord_cartesian(ylim = geo$ylims) +
    scale_fill_manual(values = cat_colors, labels = cat_labels, drop = FALSE, name = NULL) +
    theme_minimal(base_size = AXIS_TEXT_SIZE) +
    theme(
      axis.text.y = element_text(size = PATHWAY_LABEL_SIZE, lineheight = 0.95),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE),
      axis.title.x = element_text(size = AXIS_TITLE_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)
    ) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL)

  # Same axis-text/title sizes as p_bars (invisible text) so the two panels'
  # actual data areas -- not just their nominal ylim -- line up row-for-row
  # under patchwork; legend is pooled via guides="collect" below instead of
  # reserving space only inside p_bars.
  p_labels <- ggplot(geo$theme_bounds, aes(x = 0, y = y_mid)) +
    geom_text(aes(label = label), fontface = "bold", hjust = 0, size = THEME_LABEL_SIZE, lineheight = 0.95) +
    .boundary_layer(geo) +
    scale_x_continuous(limits = c(0, 1), breaks = 0.5, labels = " ") +
    # expand=c(0,0) must match p_bars' y-scale exactly (default ggplot
    # auto-expansion adds ~5% padding, which would offset this panel's rows
    # relative to p_bars' bars despite identical `ylim`).
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = geo$ylims, clip = "off") +
    .blank_axis_theme(AXIS_TEXT_SIZE, axis_text_size = AXIS_TEXT_SIZE, axis_title_size = AXIS_TITLE_SIZE)

  subtitle <- paste0(
    "Up to 8 significant pathways per biological theme (capped at 3/collection/theme), ranked by |NES|; ",
    n_themes, " of ", length(phenotype_themes), " themes matched, ", n, " pathways shown. ",
    "See scripts/utils_de_pathway.R::phenotype_themes for the exact rule."
  )
  combined <- (p_labels | p_bars) +
    plot_layout(widths = c(1.15, 4.85), guides = "collect") +
    plot_annotation(
      title = title, subtitle = paste(strwrap(subtitle, width = 110), collapse = "\n"),
      theme = theme(plot.background = element_rect(fill = "white", color = NA),
                    plot.title = element_text(size = TITLE_SIZE, face = "bold"),
                    plot.subtitle = element_text(size = SUBTITLE_SIZE))
    ) & theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "top",
              legend.text = element_text(size = LEGEND_TEXT_SIZE))

  # Dynamic dimensions: a reasonable minimum, additional height per pathway
  # row, and extra allowance per theme-group gap (both already baked into
  # `y_span` via .assign_theme_rows()'s `gap`) -- scaled up for the larger
  # fonts above so rows have comfortable, not just non-overlapping, spacing.
  y_span <- diff(geo$ylims)
  height <- max(10, 4 + y_span * 0.62)
  # Fixed but generous width: pathway labels are pre-wrapped to <=50 chars/line
  # (wrap_label() above), so the label margin's required width doesn't grow
  # with pathway count, only with the larger fonts set above.
  width <- 19
  ggsave(out_file, combined, width = width, height = height, dpi = 150, limitsize = FALSE, bg = "white")
  invisible(combined)
}

# --- 11_cross_treatment_pathway_bubble.R: pooled pathways x ALL treatments --
# Same theme-grouped, non-faceted philosophy, adapted to a dot grid: one
# shared left theme-label column, one shared pathway-name column, then one
# dot-grid panel per analysis (full / cc-independent), all row-aligned.
save_curated_pathway_bubble_figure <- function(plot_data_curated, treatments, title, out_file) {
  pathway_meta <- plot_data_curated %>%
    distinct(collection, pathway, theme) %>%
    left_join(
      plot_data_curated %>% group_by(collection, pathway) %>%
        summarise(mean_abs_nes = mean(abs(NES_full), na.rm = TRUE), .groups = "drop"),
      by = c("collection", "pathway")
    ) %>%
    .assign_theme_rows("mean_abs_nes")
  n <- nrow(pathway_meta)
  geo <- .theme_group_geometry(pathway_meta)

  d <- plot_data_curated %>%
    left_join(pathway_meta %>% select(collection, pathway, y), by = c("collection", "pathway")) %>%
    mutate(
      pathway_label = wrap_label(paste0(mapply(clean_pathway_name, pathway, collection),
                                         "  [", collection, "]"), width = 55),
      treatment = factor(treatment, levels = treatments)
    )

  nes_limit <- max(abs(c(d$NES_full, d$NES_cc_independent)), na.rm = TRUE)

  PATHWAY_LABEL_SIZE <- 4.2    # mm, geom_text
  THEME_LABEL_SIZE   <- 5.4    # mm, geom_text
  AXIS_TEXT_SIZE     <- 13     # pt, treatment names (IL/TGF/...)
  PANEL_TITLE_SIZE   <- 14     # pt, "Full analysis" / "Cell-cycle-independent analysis"
  LEGEND_TEXT_SIZE   <- 12     # pt
  LEGEND_TITLE_SIZE  <- 12.5   # pt
  TITLE_SIZE         <- 18     # pt
  SUBTITLE_SIZE      <- 12.5   # pt

  # Dot size encodes significance strength (-log10(padj), the SAME adjusted
  # p-value/FDR used everywhere else -- never raw p-value, and `padj` itself
  # is never altered). This is a DISPLAY-ONLY transform of that value for
  # legibility: the real distribution here is heavily right-skewed (checked
  # empirically -- 90% of all plotted full+cc-independent padj values across
  # both significant and non-significant dots have -log10(padj) <= ~5.3,
  # while a handful of extreme outliers reach >20), so a plain linear 0-15
  # size scale squeezes almost every real point, including most significant
  # ones, into a narrow low-size band near the non-significant minimum. Two
  # display-only fixes, applied only to this plotting variable:
  #  1. cap at 10 (data-informed: only ~4% of points exceed it) so rare
  #     extreme outliers stop stealing visual range from everything else;
  #  2. map through sqrt (`scale_size_continuous(trans="sqrt")` below) so the
  #     moderate/typical range spreads out instead of bunching near zero.
  SIG_SIZE_CAP <- 10
  make_panel <- function(nes_col, padj_col, panel_title) {
    dd <- d
    dd$NES <- dd[[nes_col]]
    dd$padj <- dd[[padj_col]]
    dd$neglog10padj <- pmin(-log10(pmax(dd$padj, 1e-20)), SIG_SIZE_CAP)
    ggplot(dd, aes(x = treatment, y = y)) +
      .band_layer(geo) +
      geom_point(aes(size = neglog10padj, color = NES)) +
      .boundary_layer(geo) +
      scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0,
                             limits = c(-nes_limit, nes_limit), name = "NES") +
      scale_size_continuous(
        trans = "sqrt", range = c(1.8, 15), limits = c(0, SIG_SIZE_CAP),
        breaks = c(0, -log10(0.05), 3, 6, SIG_SIZE_CAP),
        labels = c("NS (padj ≥ 0.05)", "0.05", "1e-3", "1e-6", "≤ 1e-10"),
        name = "padj (FDR)"
      ) +
      scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
      coord_cartesian(ylim = geo$ylims) +
      labs(title = panel_title, x = NULL, y = NULL) +
      theme_minimal(base_size = AXIS_TEXT_SIZE) +
      theme(
        axis.text.x = element_text(size = AXIS_TEXT_SIZE, face = "bold"),
        plot.title = element_text(size = PANEL_TITLE_SIZE),
        legend.text = element_text(size = LEGEND_TEXT_SIZE),
        legend.title = element_text(size = LEGEND_TITLE_SIZE),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  }
  p_full <- make_panel("NES_full", "padj_full", "Full analysis")
  p_cc   <- make_panel("NES_cc_independent", "padj_cc_independent", "Cell-cycle-independent analysis")

  # Same axis-text/title sizes as p_full/p_cc (invisible text) so all four
  # panels' actual data areas -- not just their nominal ylim -- line up
  # row-for-row under patchwork despite p_full/p_cc having a bottom legend +
  # axis text that p_names/p_themes don't.
  label_df <- d %>% distinct(y, pathway_label)
  p_names <- ggplot(label_df, aes(x = 0, y = y)) +
    geom_text(aes(label = pathway_label), hjust = 0, size = PATHWAY_LABEL_SIZE, lineheight = 0.95) +
    .boundary_layer(geo) +
    scale_x_continuous(limits = c(0, 1), breaks = 0.5, labels = " ") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = geo$ylims, clip = "off") +
    ggtitle(" ") +
    .blank_axis_theme(AXIS_TEXT_SIZE, axis_text_size = AXIS_TEXT_SIZE, title_size = PANEL_TITLE_SIZE)

  p_themes <- ggplot(geo$theme_bounds, aes(x = 0, y = y_mid)) +
    geom_text(aes(label = label), fontface = "bold", hjust = 0, size = THEME_LABEL_SIZE, lineheight = 0.9) +
    .boundary_layer(geo) +
    scale_x_continuous(limits = c(0, 1), breaks = 0.5, labels = " ") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = geo$ylims, clip = "off") +
    ggtitle(" ") +
    .blank_axis_theme(AXIS_TEXT_SIZE, axis_text_size = AXIS_TEXT_SIZE, title_size = PANEL_TITLE_SIZE)

  n_themes <- length(unique(pathway_meta$theme))
  subtitle <- paste0(
    "Biologically curated by theme (scripts/utils_de_pathway.R::phenotype_themes), up to 8 significant pathways/theme ",
    "(capped at 3/collection/theme), ranked by |NES|, pooled (", n, " unique pathways across ", n_themes, " themes) across ALL ",
    length(treatments), " treatments -- small pale dots mean that pathway simply isn't significant in that treatment, not ",
    "that it wasn't tested. NES color scale is shared between panels."
  )
  combined <- (p_themes | p_names | p_full | p_cc) +
    plot_layout(widths = c(1.3, 3.6, 1.8, 1.8), guides = "collect") +
    plot_annotation(
      title = title, subtitle = paste(strwrap(subtitle, width = 150), collapse = "\n"),
      theme = theme(plot.background = element_rect(fill = "white", color = NA),
                    plot.title = element_text(size = TITLE_SIZE, face = "bold"),
                    plot.subtitle = element_text(size = SUBTITLE_SIZE))
    ) & theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "bottom",
              legend.text = element_text(size = LEGEND_TEXT_SIZE),
              legend.title = element_text(size = LEGEND_TITLE_SIZE))

  # Dynamic dimensions, scaled up for the larger fonts above: a reasonable
  # minimum, additional height per pathway row, and the theme-group gaps
  # already baked into `y_span` via .assign_theme_rows()'s `gap`.
  y_span <- diff(geo$ylims)
  height <- max(10, 4 + y_span * 0.52)
  width <- 24
  ggsave(out_file, combined, width = width, height = height, dpi = 150, limitsize = FALSE, bg = "white")
  invisible(combined)
}
