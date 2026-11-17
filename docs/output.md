# SPADE-nf: Output

`<prefix>` is `--prefix` (default `SPADE`; CAB internal: the SRM name or `mergedSRMs_<n>`).

```
results/
├── input/            validation summary and the comparisons used
├── de/               differential expression, ERCC QC, rnk files, results object
│   └── QC/           library size, RLE, PCA (+ coordinates), ERCC QC bundle
├── enrichment/
│   ├── fgsea/<comparison>/
│   └── enrichr/<comparison>/
├── report/           HTML report, DEG summary tables, fgsea workbooks
├── merged_srms/      (CAB internal, several SRMs) merged count / TPM / FPKM / report
└── pipeline_info/    Nextflow reports, software versions
```

## input/

| File | |
|---|---|
| `<prefix>_comparisons.lst` | comparisons analysed (given, or all pairs of groups) |
| `<prefix>_validation.tsv` | samples matched, groups, species requested/detected, ERCC rows, notes |

## de/

| File | |
|---|---|
| `<prefix>_<comparison>_log2CPM.tsv` | DE table: geneID, geneSymbol, log2FC, AveExpr, t.Statistic, P.Value, FDR, group means, per-sample log2CPM |
| `<prefix>_<comparison>_t.rnk`, `_sign.log2FCxNegLog10Pval.rnk` | ranked lists for GSEA Preranked |
| `<prefix>_all_GSEAinput_TMM.txt`, `.cls` | expression and phenotype files for standard GSEA (>= 2 groups with >= 3 samples) |
| `<prefix>_SPADE_Results.rds` | R object: `QC`, `data`, `log2CPM`, `comps`, `statDEGs`, `numDEGs`, `topnDEGs`, `allDE`, `annotation`, `settings` |
| `<prefix>_ERCC_QC_stat.txt` / `.pdf` | per-sample ERCC level, normalization factors, RNA_index, ERCC_COR (ERCC mode) |
| `<prefix>_ERCC_QC_ERCC_log2count_vs_log2conc_*.png` | ERCC dose response (ERCC mode) |
| `<prefix>_MappingReport.tsv` | mapping QC (CAB internal, with a mapping report) |
| `<prefix>_log.txt` | run log, including the normalization used and why |
| `QC/<prefix>_QC_libSize.*`, `_QC_RLE*`, `_QC_PCA_TMM*` | library size, RLE and PCA plots and tables |
| `QC/<prefix>_ERCC_QC_*` | ERCC spike-in boxplots, effective library sizes, PCA (TMM vs ERCC), stat/CV/shift tables (ERCC mode) |

## enrichment/

| File | |
|---|---|
| `fgsea/<comparison>/<comparison>_<gene_set>_fgsea.tsv` | pathways with P < `--fgsea_pval_cutoff`: pval, padj, NES, size, leadingEdge, gene_ratio |
| `fgsea/<comparison>/*_fgsea_bar.pdf`, `_dot_up.pdf`, `_dot_down.pdf` | top pathways (SPADE plots) |
| `enrichr/<comparison>/<comparison>_enrichr.tsv` | Enrichr terms per library and direction |
| `enrichr/<comparison>/enrichR/*_UP_enrichR.xlsx`, `*_DOWN_enrichR.xlsx` | same, one sheet per library |
| `enrichr/<comparison>/*_{bar,dot}_{up,down}.pdf`, `*_enrichr_input_*.txt` | plots and submitted genes |

## report/

| File | |
|---|---|
| `SPADE_DE_Report_<prefix>.html` | report (CAB internal: `stdRNAseqERCC_DE_Report_<prefix>.html`) |
| `<prefix>_DEG_summary.tsv`, `<prefix>_DEG_counts.tsv` | comparisons, sample sizes and DEG numbers at the report cutoffs |
| `<prefix>_<comparison>_fgsea_results.xlsx` | fgsea results, one sheet per gene set |
| `volcanoPlots/` | one PDF per comparison (`--volcano true`) |
