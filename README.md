# SPADE-nf

**ERCC spike-in aware RNA-seq differential expression and enrichment, as a Nextflow pipeline.**

SPADE-nf runs the analysis of [SPADE](https://github.com/hongjianjin/SPADE) (Spike-in Powered Analysis
of Differential Expression and Enrichment) from the command line, on a laptop, an HPC cluster or the cloud.
It takes a gene count matrix and a sample sheet and produces:

- sample QC: library size, RLE and PCA, plus ERCC spike-in QC (level, dose response, replicate CV, global shift test);
- limma-voom differential expression for each comparison, with **ERCC normalization** when the count matrix
  has ERCC spike-ins, or library-size (TMM) normalization otherwise;
- gene set enrichment with **fgsea** on MSigDB v2026.2 (human or mouse) and, optionally, **Enrichr**;
- an HTML report and GSEA-ready `.rnk` / `.cls` files.

The analysis code is SPADE's own (copied into the container from a fixed SPADE release); the command line
and report follow the CAB RNAseqERCC pipeline.

## Pipeline

```
VALIDATE_INPUTS ─► SPADE_DE ─┬────────────────────────────────────────► SPADE_REPORT (HTML)
  (fail fast)     (QC, ERCC,  ├─► FGSEA   (comparison x gene set) ──────┤
                   limma-voom)└─► ENRICHR (per comparison, optional) ───┘
```

Changing a DEG cutoff, top-N or the report header with `-resume` re-runs only the report.

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>= 24.10) and
   Docker, Singularity or Apptainer.
2. Test on SPADE's example data (mouse, 9 samples with ERCC spike-ins):

   ```bash
   nextflow run hongjianjin/SPADE-nf -profile test,docker --outdir results_test
   ```

3. Run on your data:

   ```bash
   nextflow run hongjianjin/SPADE-nf \
       -profile docker \
       --counts counts.txt \
       --meta sample_meta.txt \
       --comparisons comparisons.lst \
       --species human \
       --outdir results
   ```

On the St. Jude HPC use `-profile stjude` (LSF + Singularity, from
[nf-core/configs](https://github.com/nf-core/configs)).

## Inputs

**Count matrix** (`--counts`, tab-delimited): raw counts, one row per gene.

| geneID | geneSymbol | bioType | annotationLevel | Sample1 | Sample2 | ... |
|---|---|---|---|---|---|---|
| ENSG00000000003.14 | TSPAN6 | protein_coding | 2 | 4100 | 4151 | ... |
| ERCC-00130 | ERCC-00130 | ERCC | 1 | 2310 | 2507 | ... |

`geneID` or `geneSymbol` is required; `bioType` and `annotationLevel` enable the protein-coding filter.
Rows whose symbol starts with `ERCC-` are spike-ins.

**Sample metadata** (`--meta`, tab-delimited):

| ID | NEWNAME | GROUP | BATCH |
|---|---|---|---|
| Sample1 | WT_1 | WT | B1 |
| Sample2 | KO_1 | KO | B1 |

`ID` must match the count-matrix columns. Optional: `COLOR`, `BATCH` or `GENDER` (added to the limma design
when they have more than one level), `SHAPE`.

**Comparisons** (`--comparisons`, optional): one `Exp_vs_Ctl` per line, e.g. `KO_vs_WT`
(log2FC = KO relative to WT). Without it, all pairs of groups are compared.

**Species** (`--species`): `human`, `mouse` or `other` (then give `--custom_gmt` for fgsea, or `--run_fgsea false`).

## Normalization

| ERCC- rows in counts | `--ercc` | Normalization |
|---|---|---|
| yes (enough reads) | `true` (default) | ERCC spike-in |
| yes | `false` | library size (TMM) |
| no, or too few ERCC reads | any | library size (TMM); the log and report say why |

## Main options

| Option | Default | |
|---|---|---|
| `--ercc` | `true` | ERCC normalization when spike-ins are present |
| `--filter_noncoding` | `true` | keep protein-coding genes only |
| `--cutoff_count`, `--cutoff_cpm` | `0`, `0` | low-expression filter; 0 = automatic (CPM 0.1) |
| `--cutoff_log2fc`, `--cutoff_pvalue`, `--cutoff_fdr` | `1`, `0.05`, `0.05` | DEG definition (report only) |
| `--topn_deg` | `10` | top up/down DEGs labelled in plots |
| `--run_fgsea`, `--fgsea_gene_sets` | `true`, `default` | MSigDB v2026.2: GO MF/BP/CC, cell type, canonical pathways, CGP, Hallmark, oncogenic (human) or tumor phenotype (mouse) |
| `--msigdb_dir` | | local MSigDB GMT files when compute nodes have no internet |
| `--run_enrichr` | `false` | Enrichr web API (needs internet) |
| `--prefix` | `SPADE` | output file prefix |

Full list: `nextflow run hongjianjin/SPADE-nf --help --show_hidden`, and [docs/usage.md](docs/usage.md).
Outputs: [docs/output.md](docs/output.md).

## Citation

If you use SPADE-nf, please cite SPADE: https://github.com/hongjianjin/SPADE, and the tools it uses:
limma/voom, edgeR, fgsea, MSigDB and Enrichr (see the report's reference list).

## License

Apache-2.0, as SPADE. MSigDB gene sets are downloaded at run time and are subject to the
[MSigDB license](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).

## Contact

Hongjian Jin, St. Jude Children's Research Hospital ([@hongjianjin](https://github.com/hongjianjin)).
