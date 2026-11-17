# SPADE-nf: Changelog

## v3.0.1 (unreleased)

First version, wrapping SPADE 3.0.1 (https://github.com/hongjianjin/SPADE).

- Command-line scripts in `bin/` (`spade_validate.R`, `spade_de.R`, `spade_fgsea.R`, `spade_enrichr.R`,
  `spade_report.R`, CAB internal `spade_merge_srms.R`) that call SPADE's analysis functions; options, outputs
  and report follow RNAseqERCC v1.1. Helpers SPADE lacks are ported from RNAseqERCC (`bin/spade_cli_helpers.R`).
- ERCC normalization when ERCC- rows are present (`--ercc`), otherwise library-size (TMM) normalization with a note.
- fgsea on MSigDB v2026.2 human/mouse GMT files (default: GO MF/BP/CC, cell type, canonical pathways, CGP,
  Hallmark, oncogenic signatures (human) / tumor phenotype (mouse)); optional Enrichr.
- SPADE-branded HTML report; CAB branding, AutoMapper inputs, multi-SRM merge and mapping QC with `--cab_internal`.
- Container `ghcr.io/hongjianjin/spade:3.0.1` (Bioconductor 3.18, R 4.3).
- Checked against RNAseqERCC v1.1 on SPADE example 3 (ERCC) and example 1 (TMM): identical DE tables.
