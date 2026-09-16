# Shared constants for the DE / pathway-analysis pipeline (scripts 05-11, 20-22).
# Source this near the top of each script (after `proj_dir` is set) instead of
# redeclaring these values locally, so a threshold change has to be made in
# only one place and can no longer silently drift between scripts.
#
# `LOG2FC_CUTOFF` is defined on DESeq2's log2FoldChange scale (base-2 log fold
# change of pseudobulk counts, script 05). limma's logFC (script 06) is on a
# different scale -- natural-log residuals after cell-cycle regression on
# Seurat's log1p-normalized `data` layer -- so applying "the same" minimal-
# effect gate to limma results requires a scale conversion. Use
# `fc_cutoff("limma")` rather than `LOG2FC_CUTOFF` directly when filtering
# limma output; use `fc_cutoff("deseq2")` (equal to `LOG2FC_CUTOFF`) for
# DESeq2 output, so the conversion lives in exactly one place.

SIG_ALPHA <- 0.05
LOG2FC_CUTOFF <- 0.1
PATHWAY_TOP_N <- 15
# Size of the experiment-derived "induced" reference signatures (scripts 20/21:
# top N genes ranked by effect size among those passing SIG_ALPHA + fc_cutoff()
# + positive direction). Deliberately separate from PATHWAY_TOP_N -- these are
# different kinds of top-N lists (genes for a reference signature vs. pathways
# for a plot) and must not silently share one value.
EXPERIMENT_SIGNATURE_TOP_N <- 50
CONDITION_LEVELS <- c("Control", "IL", "TGF", "PDGF_B", "PDGF_C")
TREATMENTS <- setdiff(CONDITION_LEVELS, "Control")
REPLICATE_LEVELS <- c("r1", "r2", "r3")
COLLECTIONS <- c("Hallmark", "GO_BP", "Oncogenic", "Reactome", "KEGG")

fc_cutoff <- function(engine = c("deseq2", "limma")) {
  engine <- match.arg(engine)
  if (engine == "deseq2") LOG2FC_CUTOFF else LOG2FC_CUTOFF * log(2)
}
