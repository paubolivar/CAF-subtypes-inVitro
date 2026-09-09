# Project Context

## Experimental system

In vitro cancer-associated fibroblasts (CAFs), CAF2 cell line.

## Cell type

CAF2 cells only. Other cell types are not expected; unexpected populations
should initially be treated as possible contamination, not assumed to be
real biological populations.

## Technology / platform

10x Genomics Flex (Fixed RNA Profiling), processed with Cell Ranger multi
v10.0.0. Gene Expression only (no VDJ, antibody capture, or CRISPR).
Samples are demultiplexed by probe barcode. Reference: GRCh38-2024-A.

## Sample structure

15 samples total, all pooled into a single GEM well / sequencing run
(one technical batch).

## Biological replicates

3 independent biological replicates (r1, r2, r3) per treatment condition.

## Treatment conditions

Control, IL, TGF, PDGF_B, PDGF_C (5 conditions).

## Biological question

How the treatments alter CAF transcriptional states.

## Other important experimental facts

- A second cell line ("Mel") is planned per the lab's sample-allocation
  spreadsheet but has no sequencing data in `data/` yet.
- The specific interleukin behind the "IL" condition label is not yet
  confirmed.
