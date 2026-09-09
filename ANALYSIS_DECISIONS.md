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
