# Analysis Decisions

## 2026-09-09 — Statistical unit for treatment comparisons

### Decision
The biological replicate, not the individual cell, is the unit of
inference for treatment-level statistical comparisons.

### Rationale
r1, r2, and r3 are confirmed independent biological replicates (separate
wells/passages, independently treated). Cells within a replicate are not
independent of each other, so treating individual cells as replicates
would inflate the effective sample size.

## 2026-09-09 — Analysis scope

### Decision
Analyze the CAF2 line now. The Mel cell line will be handled as a
separate future addition once its sequencing data is available.

### Rationale
Only CAF2 sequencing data currently exists in `data/`; Mel is defined in
the planning spreadsheet but has not been sequenced/added yet.

## 2026-09-09 — Analysis stack

### Decision
Use R / Seurat for downstream analysis.

## 2026-09-10 — QC filtering applied for treatment-response analysis

### Decision
For the global treatment-response pseudobulk analysis, actually apply the 5-MAD
per-sample `nCount_RNA` outlier filter and `percent.mt > 5` filter that
`01_qc_explore.R` had only computed descriptively. Clusters 15/16 (flagged in
`03_qc_deepdive.R` as possibly low-quality/artifactual) are left untouched, and no
automated classifier (e.g. scDblFinder) was run.

### Rationale
These thresholds were already computed and reviewed; applying them is a light,
well-understood filter (0.47% of cells removed overall, evenly across samples).
Clusters 15/16 and doublet calling remain open questions not yet resolved with the
user, so they were deliberately excluded from this filtering step.

### Consequence
A new object, `results/caf2_merged_filtered_ccscored.rds`, was created
(`scripts/04_build_qc_filtered_object.R`) distinct from the unfiltered objects saved
by `01_`/`02_`. All treatment-response results (Parts A/B/C) are built on this
filtered object. If clusters 15/16 are later confirmed as artifacts, this analysis
would need to be re-run after excluding them.

## 2026-09-10 — Cell-cycle-independent method for treatment response

### Decision
For the cell-cycle-independent treatment-response analysis (Part B), regress
`S.Score`/`G2M.Score` out of each cell's log-normalized expression via a vectorized
OLS projection (mathematically equivalent to Seurat's `ScaleData(vars.to.regress=...,
do.scale=FALSE, do.center=FALSE)` linear-model residuals), then mean-aggregate the
residuals per sample into a pseudobulk matrix and test with `limma` rather than
DESeq2.

### Rationale
The tool choice for each part follows directly from what kind of data that part's
pseudobulk matrix actually contains — it was not an independent, arbitrary choice
between DESeq2 and limma.

**Part A (full, cell cycle included)**: the pseudobulk input is raw summed UMI
counts per sample (`AggregateExpression` summing counts across cells within each of
the 15 samples). These are still non-negative integers, so DESeq2's negative-binomial
GLM — built specifically to model count data and its characteristic mean-variance
relationship — is the standard, appropriate tool. This is the conventional pseudobulk
DGE approach.

**Part B (cell-cycle-independent)**: to remove the cell-cycle contribution,
`S.Score`/`G2M.Score` are regressed out of each cell's *log-normalized* expression
value (fitting `expression ~ S.Score + G2M.Score` per gene and keeping the
residuals). Residuals from a linear regression are continuous numbers — they can be
negative, they are not counts, and they no longer follow a count distribution. Once
those residuals are mean-aggregated per sample into a pseudobulk matrix, the result
is a continuous-valued matrix, not a count matrix. DESeq2 requires actual counts (it
would error, or produce meaningless results, on negative/non-integer input) — so it
is not usable here. `limma` (`lmFit`/`eBayes`), which fits linear models to
continuous, roughly-Gaussian data, is the standard tool for that kind of matrix (it
is also the conventional choice for microarray/continuous expression data more
generally).

In short: Part A stayed on counts throughout, so it kept a count-based test
(DESeq2). Part B's cell-cycle-removal step converts the data from counts to
continuous residuals *before* pseudobulk aggregation, so it needed a test built for
continuous data (limma). The two parts therefore use different differential-testing
frameworks not because of a preference between DESeq2 and limma in the abstract, but
because the chosen cell-cycle-removal method (regress at the single-cell level, then
aggregate) changes the data type partway through the Part B pipeline only.

(Separately, Seurat's per-gene `ScaleData` regression was impractically slow at ~96k
cells — over 3h single-threaded, and `future` parallelization failed/was inefficient
due to whole-object serialization to workers — so a mathematically-equivalent
vectorized matrix projection was used instead to compute the same residuals faster.
This is an implementation detail and does not change the statistical framework
described above.)

### Consequence
Part A (DESeq2/counts) and Part B (limma/continuous residuals) are not on identical
statistical footing — different null models, variance-shrinkage, and power. Part C's
"cell-cycle-driven vs. robust" classifications therefore partly reflect this
methodological gap, not purely cell-cycle removal; this is stated explicitly in
`results/pseudobulk_cc_independent/interpretation_cc_independent.md` and
`results/comparison_full_vs_cc_independent/interpretation_comparison.md`. A recommended follow-up is a
DESeq2-covariate sensitivity re-run (`~S.Score+G2M.Score+condition`) to isolate the
two effects.

## 2026-09-10 — FGSEA gene set collections

### Decision
Use `msigdbr` (Homo sapiens) with five MSigDB collections for pathway analysis:
Hallmark (H), GO Biological Process (C5/GO:BP), Oncogenic Signatures (C6), Reactome
(C2/CP:REACTOME), and KEGG (C2/CP:KEGG_LEGACY).

### Rationale
These give broad, standard coverage (curated hallmark processes, GO ontology,
oncogenic programs, and two canonical-pathway resources) for a first pathway-level
look, since no pathway analysis existed in this project before.

### Consequence
The installed `msigdbr` version (25.1.1) splits KEGG into `CP:KEGG_LEGACY` and
`CP:KEGG_MEDICUS`; `CP:KEGG_LEGACY` was used as the closer analog to classic KEGG
pathways. Re-running with a different `msigdbr` version should re-check this label
via `msigdbr_collections()`.

## 2026-09-10 — Treatment contrasts for global treatment-response analysis

### Decision
Exactly 4 contrasts, each comparing one treatment to Control: Control vs IL, Control
vs TGF, Control vs PDGF_B, Control vs PDGF_C. Control is always the reference level
in both DESeq2 and limma design matrices.

### Rationale
Matches the experimental design (one control, four treatments) and the biological
question (how each treatment alters CAF2 transcriptional state relative to Control).

### Consequence
All Part A/B/C outputs are organized per-contrast (`Control_vs_<TRT>`), giving a
fixed 4 contrasts x 5 collections x 2 parts output matrix, useful for checking
completeness of results.

## 2026-09-11 — Corrected significance testing for published CAF signature scoring

### Decision
Significance testing for the published-signature scoring analysis
(`results/caf_signature_scoring/`, `results/caf_signature_scoring_zscore/`) was
changed from a per-cell Wilcoxon rank-sum test to `limma` on per-replicate means
(15 pseudobulk points, n=3/group vs Control), matching the pseudobulk DGE approach
used elsewhere in this project.

### Rationale
The original per-cell test treated each cell as an independent sample, which
conflicts with the 2026-09-09 decision above (biological replicate is the unit of
inference) — cells from the same replicate/well are correlated, not independent.
This was a real pseudoreplication problem, not just "large N inflates
significance": with tens of thousands of correlated cells per group, nearly every
test returned p≈0 regardless of effect size, making the significance values
uninterpretable.

### Consequence
`signature_treatment_summary.csv` in both output folders now contains
replicate-level `padj` (limma, BH-adjusted). The original per-cell version is
preserved for reference as `signature_treatment_summary_percell_exploratory.csv`
in each folder but should not be used for significance claims. Effect-size (delta)
values were materially unchanged between the two versions, so prior biological
interpretations based on delta remain valid; only the significance values needed
correcting.

## 2026-09-11 — Fixed inconsistent UMAP embeddings between population-structure and signature-scoring figures

### Decision
Clustering/UMAP for reporting purposes is now computed directly on the
QC-filtered, cell-cycle-scored object (`scripts/23_cluster_qc_filtered.R`),
not on the unfiltered object.

### Rationale
`02_cluster_explore.R` clustered the *unfiltered* merged object (96,201
cells) — that was correct for its original purpose, an exploratory
pre-filtering QC step used to flag clusters 15/16. `14_compute_umap_qc_filtered.R`
separately computed a UMAP on the *filtered* object (95,750 cells) for
signature-scoring feature plots. Because these were fit on different cell
populations, they are not the same embedding — the population-structure
figure and the signature-scoring UMAPs in the treatment-response report
visibly disagreed. This was caught by manual review, not by any automated
check.

### Consequence
`23_cluster_qc_filtered.R` reclusters the filtered object (same
FindVariableFeatures -> ScaleData -> RunPCA -> RunUMAP recipe and seed as
before, plus FindNeighbors/FindClusters/FindAllMarkers), overwriting
`caf2_merged_filtered_ccscored_umap.rds` in place. `02_cluster_explore.R`
and `results/clustering/` are left untouched as the historical record of the
pre-filtering QC exploration; `14_compute_umap_qc_filtered.R` is marked
superseded. Cluster ID numbers shifted slightly on re-clustering (e.g. the
IL-dominant islands are now clusters 6/10/13, not 6/10/12) — the underlying
biology and marker genes are unchanged, only the labels. All report files
(then `results/report/`, now `results/CAF2/report/` after the reorganization
below) were corrected accordingly. While re-verifying this fix, an unrelated
reporting error was also caught and corrected: the report had stated "~89,750
of ~90,200 cells retained" after QC filtering; the actual `qc_filter_audit.csv`
totals are 95,750 of 96,201 (0.47% removed).

## 2026-09-11 — Multi-cell-line pipeline generalization

### Decision
The Mel cell line's sequencing data arrived (same experimental design as CAF2:
5 conditions x 3 replicates, same 10x Flex platform, same gene panel — verified
by direct comparison). Rather than duplicating scripts per cell line, every
pipeline script (`01`, `04`-`12`, `15`-`23`) was made to accept an optional
`CELL_LINE` command-line argument (default `CAF2`, preserving prior behavior),
and all `results/` output paths were changed from flat
(`results/<analysis_name>/`) to per-cell-line
(`results/<CELL_LINE>/<analysis_name>/`). CAF2's existing outputs were moved
(not recomputed) into `results/CAF2/` to match. A new
`scripts/run_pipeline.sh <CELL_LINE>` runs the full ordered pipeline in one
call; see `scripts/README_pipeline.md` for the exact step order and
per-cell-line usage instructions.

### Rationale
The project's own experimental design (`PROJECT_CONTEXT.md`) always
anticipated more than one cell line. Hardcoding "CAF2" into output paths and
sample-name parsing would have meant either duplicating ~19 scripts per cell
line (error-prone — any future fix would need to be applied N times) or
silently overwriting CAF2's results when analyzing Mel. Parameterizing once,
now, keeps a single source of truth for the analysis logic while giving each
cell line a fully independent results tree.

### Consequence
- Results are never pooled or numerically compared across cell lines by this
  pipeline; each cell line's `results/<CELL_LINE>/` tree is self-contained,
  including its own reference-signature panel
  (`results/<CELL_LINE>/experiment_signatures/`, source-tagged
  `<CELL_LINE>exp`/`<CELL_LINE>expFull` rather than `CAF2exp`/`CAF2expFull`).
- `02_cluster_explore.R` and `03_qc_deepdive.R` remain CAF2-specific, one-off
  historical QC-exploration scripts (not parameterized, not part of the
  standard pipeline) — their original purpose (investigating clusters 15/16 on
  the unfiltered CAF2 object) doesn't generalize to other cell lines the way
  the validated pipeline steps do; see README_pipeline.md.
- The published CAF gene-signature list
  (`data/published_data/CAF_gene_signatures/CAF_gene_signatures_combined.csv`,
  from `13_combine_caf_gene_signatures.R`) is shared, unchanged, across cell
  lines — it is external reference data, not cell-line-specific.
- Mel's pipeline run and its own interpretation/report are a direct,
  independent application of this same methodology to new data — any
  biological conclusions for Mel must come from Mel's own results, not be
  assumed from CAF2's.

## 2026-09-11 — Cross-cell-line comparison scope and layer choice

### Decision
A dedicated cross-cell-line comparison (`scripts/compare_cell_lines.R`,
`results/compare_CAF2_vs_Mel/CAF2_vs_Mel_comparison_report.Rmd`) compares CAF2
and Mel's already-computed results directly: gene/pathway overlap (intersection,
Jaccard) and signature-score correlation. It performs no new statistical tests
beyond these overlap/correlation summaries, and uses the **cell-cycle-independent
layer only** for every gene- and pathway-level comparison, not the full
(cell-cycle-included) layer.

### Rationale
Both cell lines' cell-cycle-independent values come from the identical
limma-on-CC-regressed-residuals method on the same expression scale, making them
directly comparable in a way the full-analysis DESeq2 values (fit per cell line
on that cell line's own raw counts) are not guaranteed to be. Showing both layers
side by side for two cell lines would also have produced a 4-panel comparison per
figure, judged too dense to be useful.

### Consequence
`results/compare_CAF2_vs_Mel/` findings: population structure (IL/TGF islands
vs. PDGF non-separation) and TGF's myofibroblast gene program (2,821 shared
significant cell-cycle-independent genes, Jaccard 0.376) replicate strongly
between cell lines; IL also overlaps strongly (Jaccard 0.395). PDGF_C is the
clear exception (Jaccard 0.02, only 90 shared genes) — driven mainly by CAF2
having a much smaller significant-gene pool for PDGF_C (139 genes) than Mel
(4,559 genes), not by the shared genes disagreeing in direction. PDGF_B overlaps
moderately (Jaccard 0.271, 2,129 shared genes) despite strong signature-level
correlation (r = 0.90, 0.87 for PDGF_B/PDGF_C) — read together, this suggests
PDGF's response program is broadly conserved but PDGF_C specifically is a much
weaker, narrower effect in CAF2 than in Mel, elaborated as an explicitly-labeled
hypothesis (not a firm conclusion) in the comparison report. Also notable: Mel's
PDGF-induced proliferation/mitotic pathways remain significant in the
cell-cycle-independent layer where CAF2's do not — two competing explanations
(incomplete cell-cycle correction vs. a genuinely stronger proliferative
response in Mel) are both stated as unresolved by this data.
*(Numbers in this entry were updated 2026-09-12 after the `robust_independent`
gene classification below was removed — see that entry. The overlap numbers
above reflect the post-removal significance-based gene comparison; the
qualitative conclusions were unchanged, but PDGF_B's Jaccard moved from 0.057 to
0.271 and PDGF_C's from 0.010 to 0.02 once genes excluded only for failing the
old 50%-effect-size-retention rule were restored to the comparison.)*

## 2026-09-12 — Removed the `robust_independent` gene classification (50%-effect-size-retention rule)

### Decision
Removed the gene-level classification previously computed in
`scripts/07_compare_full_vs_cc_independent.R`, which required a gene to be (a)
significant in both the full and cell-cycle-independent analyses (padj < 0.05
each), (b) concordant in direction, and (c) retain at least 50% of its full-
analysis effect size in the cell-cycle-independent analysis
(`ATTENUATION_RATIO = 0.5`) before being labeled `robust_independent` and used
downstream. Genes are now selected and ranked purely by significance in the
cell-cycle-independent analysis (`adj.P.Val_cc_independent < 0.05`) and
`|logFC_cc_independent|` — no comparison to the full analysis, no effect-size-
retention ratio. This selection logic is now applied consistently everywhere a
"top genes" list is produced: the cross-treatment DGE heatmap (script 10), the
experiment-derived reference signatures (scripts 20/21), the per-treatment
significant-gene exports (script 22), and the cross-cell-line gene comparison
(`compare_cell_lines.R`). Pathway-level classification
(`sig_both_concordant`/`sig_full_only`/`sig_cc_independent_only`/
`not_sig_either`) is unaffected — it was always based on significance and NES
direction only, with no effect-size ratio, so it needed no change.

### Rationale
The user flagged this threshold as biologically arbitrary during review: a gene
like `IL6`, whose log2FC dropped from ~5.5 in the full analysis to ~2.35 in the
cell-cycle-independent one, is still a real ~5-fold change by normal RNA-seq
standards — not a non-finding — yet the 50% ratio excluded it from the
"robust" gene set and therefore from every downstream "top genes" list,
heatmap, and reference signature built from it. The round-number 50% cutoff had
no mathematical or biological grounding (discussed conceptually with the user
before any change was made); it was carried over informally from the earlier
plan without a principled basis. Significance plus effect-size ranking, with no
retention-ratio gate, is a simpler and more defensible standard that keeps
genes with large, real, still-significant effects instead of discarding them
for an arbitrary reason.

### Consequence
- `gene_concordance_Control_vs_<TRT>.csv` (script 07 output) no longer has a
  `category` column; it is a plain joined table of full- and cell-cycle-
  independent-analysis statistics per gene. `gene_category_counts.png` was
  removed (no longer meaningful without gene categories).
- Several conclusions drawn earlier in this project (before this removal, when
  writing the CAF2-vs-Mel comparison) turned out to be artifacts of the 50%
  threshold rather than real biology, and were corrected in place rather than
  silently changed: Mel's IL cytokine genes (`CXCL8`, `CXCL1`, `MMP1`, `MMP3`,
  etc.) and PDGF_B ECM-suppression genes were previously described as failing
  to "replicate robustly" in Mel — they were always significant in Mel's
  cell-cycle-independent analysis, just excluded from the old "robust" set for
  falling under the ratio. CAF2's PDGF_C experiment-derived signature grew from
  4 genes to the full 20 once the ratio gate was removed. The genuine,
  still-persisting cross-cell-line difference is PDGF_C's much smaller overall
  significant-gene pool in CAF2 (139 genes) vs. Mel (4,559 genes) — a real
  quantitative difference, not a classification artifact.
- Every downstream document was updated to match: both cell-line Rmd/`.md`
  reports (`results/CAF2/report/`, `results/Mel/report/`), both
  `interpretation_comparison.md` files, both `experiment_signatures/README.md`
  files, the cross-cell-line comparison report and its
  `ANALYSIS_DECISIONS.md` entry above (2026-09-11), and the published CAF2
  artifact HTML (`https://claude.ai/code/artifact/22aa93c0-6097-4be4-9aa5-9d230b3b395d`,
  now Version 6).
