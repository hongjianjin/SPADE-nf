# SPADE-nf: Changelog

## v3.0.1

**Container rebuilt (2026-11-17).** The first `ghcr.io/hongjianjin/spade:3.0.1` image lacked the R package
`ggpubr`, so `SPADE_DE` failed with `there is no package called 'ggpubr'`. The image was rebuilt under the
same tag from the fixed Dockerfile on `main` (CRAN pinned to a 2024-04-24 snapshot; the build now fails if a
required package does not load). If you pulled the image before the rebuild, delete the cached copy, e.g.
`rm $NXF_SINGULARITY_CACHEDIR/ghcr.io-hongjianjin-spade-3.0.1.img`. To rebuild the image yourself, use the
Dockerfile from `main`, not from the `v3.0.1` tag.

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
