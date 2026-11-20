# Weed Rhizome dsRNA & RT-qPCR Primer Design (`Biostrings`)

An R and `Biostrings` workflow for evaluating coding sequences (CDS) from a whole-genome weed annotation, ranking candidate targets by 21-mer uniqueness, scoring sliding windows for double-stranded RNA (dsRNA) construct design, and designing flanking RT-qPCR primers outside the dsRNA target region to measure transcript knockdown.

## Repository Structure

```text
Weedgenes/
├── README.md
├── .gitignore
└── scripts/
    ├── 01_rank_cds_21mer_uniqueness.R
    ├── 02_score_dsrna_candidates.R
    └── 03_design_flanking_qpcr_primers.R

Rscript scripts/01_rank_cds_21mer_uniqueness.R
Rscript scripts/02_score_dsrna_candidates.R
Rscript scripts/03_design_flanking_qpcr_primers.R

