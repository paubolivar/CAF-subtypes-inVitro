# Project Context

## Experimental system

In vitro cancer-associated fibroblasts (CAFs). Two independent cell lines, **CAF2** and **Mel**, are analyzed in parallel with the identical experimental design.

Each cell line is its own independent biological system. Results are **not pooled across cell lines**, and a finding in one cell line is not assumed to hold in the other unless separately confirmed.

## Cell type

CAF2 and Mel cells, each analyzed separately. Other cell types are not expected within either line's data; unexpected populations should initially be treated as possible contamination, technical artifacts, or unexpected states, not assumed to be real biological populations.

## Technology / platform

10x Genomics Flex (Fixed RNA Profiling), processed with Cell Ranger multi v10.0.0.

Gene Expression only (no VDJ, antibody capture, or CRISPR).

Samples are demultiplexed by probe barcode.

Reference: GRCh38-2024-A.

Same platform and reference for both cell lines.

## Sample structure

15 samples per cell line (30 total across CAF2 and Mel), with each cell line's 15 samples pooled into its own single GEM well / sequencing run.

This represents one technical batch per cell line.

Raw per-sample Cell Ranger outputs:

`data/<CELL_LINE>_results/per_sample_outs/<CELL_LINE>_<replicate>_<condition>/`

## Biological replicates

3 independent biological replicates (r1, r2, r3) per treatment condition.

The biological replicate is the experimental unit for treatment-level inference.

## Treatment conditions

* Control
* IL
* TGF
* PDGF_B
* PDGF_C

The specific interleukin represented by the "IL" condition label is not yet confirmed for either cell line.

## Biological question

How the treatments alter CAF transcriptional states.

The primary biological variable is **treatment**.

The main analytical goal is to distinguish:

1. global transcriptional treatment responses,
2. treatment responses associated with cell-cycle activity,
3. treatment responses that remain after accounting for cell-cycle effects, and
4. heterogeneity of CAF responses within each treatment.

---

# Environment / software

All analysis should be performed inside the project's designated conda environment, currently `seurat5_scRNA`.

The environment has been confirmed to contain:

* Seurat
* DESeq2
* limma
* edgeR
* msigdbr
* fgsea
* apeglm
* ashr

Do not install packages into the system R/Python environment.

Do not create a new conda environment unless explicitly requested.

Before running analyses, verify that the correct conda environment is active.

If a required package is missing, report this before installing anything.

Prefer adding required packages to the existing project environment rather than creating a new environment.

---

# Analysis execution rules

These rules are important for maintaining a controlled and reproducible analysis.

### Follow the agreed analysis plan

Do not advance to later stages of the project unless explicitly asked.

Do not add new analyses, filtering strategies, statistical procedures, gene-selection procedures, or alternative workflows simply because they appear potentially useful.

If an additional analysis seems useful, **propose it first and wait for approval**.

### Do not silently modify the analysis

Do not introduce additional filtering of genes or cells unless that filtering criterion has been explicitly agreed upon.

In particular, do not:

* filter DE genes based on arbitrary fold-change or expression thresholds,
* remove genes because they are "uninformative" without an agreed rationale,
* restrict pathway analysis to a manually selected subset of DE genes,
* introduce additional significance thresholds beyond those specified in the analysis plan,
* remove treatments, replicates, or cells because they produce inconvenient results,
* or change the statistical model because another model produces a more interpretable result.

Any filtering or transformation that materially changes the biological result must be explicitly described and approved before implementation.

### Preserve the distinction between exploratory and primary analyses

Exploratory analyses may be useful for understanding the data, but they should not automatically become part of the primary analysis pipeline.

Clearly label exploratory analyses as exploratory.

Do not use exploratory findings to retroactively redefine the primary analysis without discussion.

### Before substantial analysis

For any new analytical stage:

1. State what will be done.
2. State which part of the project plan it corresponds to.
3. Identify important assumptions.
4. Identify any consequential analytical choices.
5. Mention a reasonable alternative if the choice could materially affect interpretation.
6. Wait for approval before implementing a new substantial analysis.

For small routine steps that directly implement an already-approved analysis, approval is not required.

### Protect the biological question

Do not optimize the analysis toward obtaining expected biological results.

Unexpected or negative results are valid results and should not be removed or modified merely because they are biologically inconvenient.

---

# Current project status

The project is currently being developed primarily using **CAF2**.

The Mel data have now arrived and should be kept as a separate cell-line analysis. Do not automatically expand the current CAF2 analysis to Mel unless explicitly requested.

For CAF2, initial exploratory QC, normalization, clustering, and cluster-level marker analyses have already been performed.

The current work is focused on establishing a defensible QC/filtering decision before proceeding with the primary treatment-response analysis.

The intended next major analysis is **pseudobulk treatment-response DGE**, followed by pathway analysis.

Do not proceed to new treatment-response analyses until the current QC decisions have been reviewed and approved.

---

# Project plan

This project plan can evolve, but changes to the plan should be explicit and documented.

## 1. QC / CAF identity

* QC and filtering
* Confirm CAF identity
* Assess sample / replicate quality
* Assess potential technical artifacts and suspicious populations
* Assess potential doublets, without automatically removing cells based on a doublet score
* Establish and document final QC decisions

## 2. Global treatment response

### A. FULL treatment-response analysis

* Pseudobulk DGE
* Biological replicate as experimental unit
* Control vs each treatment
* Cell cycle INCLUDED
* FGSEA / pathway analysis
* Pathway enrichment plots
* Brief biological interpretation

### B. CELL-CYCLE-INDEPENDENT treatment-response analysis

* Pseudobulk DGE
* Biological replicate as experimental unit
* Control vs each treatment
* Cell-cycle contribution REMOVED / accounted for using an explicitly defined model
* FGSEA / pathway analysis
* Pathway enrichment plots
* Brief biological interpretation

### C. Compare the two treatment-response analyses

* Which treatment effects are strongly associated with cell-cycle activity?
* Which treatment effects remain after accounting for cell cycle?
* Which biological responses are robust across both analyses?
* Which pathways appear primarily cell-cycle-associated?
* Which pathways become more apparent after accounting for cell cycle?

The comparison should be interpreted biologically rather than treated as a simple list of "significant in one / significant in the other" genes.

## 3. CAF heterogeneity within treatment

Analyze each treatment separately.

* Transcriptional dispersion / heterogeneity
* Distribution of continuous CAF programs
* Distribution / abundance of discrete states
* Identify treatment-specific CAF axes / states
* Determine whether treatment changes the internal CAF landscape
* Identify candidate genes associated with heterogeneity axes or discrete states

Any dimensionality reduction, clustering, trajectory, or state-inference method should be evaluated for suitability to the actual experimental design before implementation.

## 4. Integrative interpretation

* Treatment → global CAF response
* Treatment → CAF programs
* Treatment → heterogeneous CAF responses
* Treatment × CAF state/program
* Cell-cycle-dependent vs cell-cycle-independent effects
* Identify candidate mechanisms
* Prioritize candidate markers / gene signatures

---

# Important biological interpretation principles

* Treatment is the primary biological variable for treatment-response analyses.
* Biological replicates are the experimental units for treatment-level inference.
* Cell-cycle activity may represent genuine treatment biology and should not automatically be treated as technical noise.
* Cell-cycle regression/accounting should therefore be treated as a biological comparison, not automatically as a "correction."
* CAF state labels such as iCAF or myCAF should not be assigned solely from treatment identity.
* Known CAF signatures should be evaluated against the observed transcriptional states rather than used to define those states in advance.
* Unexpected transcriptional populations should be investigated before being assigned biological identities.
* A computational result should not be described as a biological mechanism unless the available data support that interpretation.

## Other important experimental facts

* The Mel cell line's sequencing data arrived 2026-09-11. The design is the same 5 conditions × 3 biological replicates.
* Analysis scripts are parameterized by cell line (see `scripts/README_pipeline.md`) so the same pipeline can be re-run identically for CAF2, Mel, or any future cell line added the same way.
* Results are kept in fully independent per-cell-line folders (`results/<CELL_LINE>/...`).
* Cell lines must never be pooled or numerically compared without an explicit, separate cross-cell-line analysis step.

## Do not invent analytical criteria

Do not introduce new numerical cutoffs, thresholds, ranking criteria, "robustness" criteria, gene-selection rules, or other decision rules that are not explicitly specified in the approved analysis plan.

This includes seemingly reasonable criteria such as retention ratios, minimum effect-size changes, numbers of genes retained, or custom definitions of "robust" results.

If a result requires such a criterion, do not choose one independently. Describe the issue, explain the possible options and their implications, and ask for approval before applying the criterion.

Do not use an exploratory or informal criterion from an earlier discussion as though it were an established methodological decision unless it is explicitly documented as an approved decision.

A gene or pathway should not be excluded from downstream analysis solely because it fails an invented or newly introduced criterion.


## Data provenance and plotting integrity — NON-NEGOTIABLE

Primary statistical result tables are the source of truth.

This includes outputs from:

* DESeq2
* limma
* FGSEA / MSigDB pathway analysis
* other primary statistical analyses in the pipeline

### Primary results must never be modified for plotting

Downstream scripts must NOT modify, overwrite, or restructure primary statistical result files simply to generate figures.

In particular:

* Do not change statistical values.
* Do not rename native DESeq2 or limma result columns.
* Do not overwrite primary result files with filtered/selected versions.
* Do not manually add genes or pathways to statistical results.
* Do not manufacture rows for plotting.
* Do not assign or alter logFC, NES, p-values, padj/FDR, or other statistical values manually.

Filtering, ranking, labeling, annotation, biological curation, and top-N selection must occur downstream and affect only what is DISPLAYED.

**Selection applies to plots, not to the underlying statistical results.**

### Every plotted element must have data provenance

Every gene or pathway appearing in a plot must be traceable to an actual row in the appropriate statistical result.

For gene-level figures:
`DESeq2/limma result → significance/filtering → plot selection → plotted gene`

For pathway figures:
`FGSEA result → FDR significance → optional biological selection → abs(NES) ranking → plotted pathway`

No gene or pathway may appear simply because it was manually specified as biologically interesting.

### Curated pathway plots

Biological curation is allowed ONLY as a selection mechanism.

Theme definitions and keyword/regex catalogs may identify relevant pathways, but they must only select pathways that actually exist in the FGSEA results.

Therefore:

**Biological curation may select from the data; it may never manufacture the data.**

A curated pathway must:

1. exist in the corresponding FGSEA/MSigDB result;
2. satisfy the required FDR significance criterion;
3. retain its original NES and statistical values;
4. be selected/ranked according to the documented plotting rules.

If no significant pathway exists for a biological theme, nothing should be displayed for that theme. Do not force representation of a biological theme.

Names such as `Interferon_Response`, `TNF_NFkB`, `ECM_Organization`, etc. are our biological annotation/grouping categories. They must never be presented as if they were original MSigDB pathway names.

### FULL vs CC-independent

When FULL and CC-independent results are combined for visualization, provenance must remain explicit.

For every displayed value, it must be possible to determine:

* whether it originated from FULL or CC-independent;
* its original result row;
* original NES/logFC;
* original padj/FDR;
* significance category;
* why that value was selected for display.

Never silently substitute values between FULL and CC-independent analyses.

### General rule when modifying plotting code

Whenever modifying or creating a figure, verify that:

**result table → filtering → selection → plotting dataframe → figure**

is fully traceable.

Manual lists may be used for filtering, biological grouping, annotation, or labeling, but never to create statistical observations that are absent from the source results.

Existing DESeq2, limma, FGSEA/MSigDB, and other primary statistical outputs must remain untouched unless I explicitly request a change to the statistical analysis itself.

## Project plan

This project pan can evolve.

CAF TREATMENT PROJECT
│
├── 1. QC / CAF identity
│   ├── QC and filtering
│   ├── Confirm CAF identity
│   └── Assess sample / replicate quality
│
├── 2. Global treatment response
│   │
│   ├── A. FULL treatment-response analysis
│   │   ├── Pseudobulk DGE
│   │   ├── Control vs each treatment
│   │   ├── Cell cycle INCLUDED
│   │   ├── FGSEA / pathway analysis
│   │   ├── Pathway enrichment plots
│   │   └── Brief biological interpretation
│   │
│   ├── B. CELL-CYCLE-INDEPENDENT treatment-response analysis
│   │   ├── Pseudobulk DGE
│   │   ├── Control vs each treatment
│   │   ├── Cell-cycle contribution REMOVED / accounted for
│   │   ├── FGSEA / pathway analysis
│   │   ├── Pathway enrichment plots
│   │   └── Brief biological interpretation
│   │
│   └── C. Compare the two treatment-response analyses
│       ├── Which treatment effects are driven by cell cycle?
│       ├── Which treatment effects remain independent of cell cycle?
│       └── Identify robust biological treatment responses
│

├── 4. CAF heterogeneity WITHIN treatment
│   │
│   ├── Analyze each treatment separately
│   ├── Transcriptional dispersion / heterogeneity
│   ├── Distribution of continuous CAF programs
│   ├── Distribution / abundance of discrete states
│   ├── Identify treatment-specific CAF axes / states
│   └── Determine whether treatment changes the CAF
│       landscape internally
│
└── 5. Integrative interpretation
    ├── Treatment → global CAF response
    ├── Treatment → CAF programs
    ├── Treatment → heterogeneous CAF responses
    ├── Treatment × CAF state/program
    ├── Cell-cycle-dependent vs cell-cycle-independent effects
    ├── Identify candidate mechanisms
    └── Prioritize candidate markers / gene signatures