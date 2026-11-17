# SPADE-nf: Usage

## Running the pipeline

```bash
nextflow run hongjianjin/SPADE-nf -r v3.0.1 \
    -profile <docker|singularity|apptainer|stjude> \
    --counts counts.txt --meta sample_meta.txt [--comparisons comparisons.lst] \
    --species <human|mouse|other> --outdir results
```

Use `-r <version>` for reproducible runs and `-resume` to reuse finished steps. Options can also go in a
YAML file passed with `-params-file params.yaml`.

## Inputs

See the [README](../README.md#inputs) for the count matrix, metadata and comparison formats. Checks done by
`VALIDATE_INPUTS` before any analysis:

- metadata has `ID`, `NEWNAME`, `GROUP`; IDs and names are unique and use only `A-Z a-z 0-9 . _ -`;
- metadata IDs match count-matrix columns (unmatched IDs are reported and dropped);
- `--species` agrees with the gene IDs (`ENSG` / `ENSMUSG`) or the case of the gene symbols;
- every comparison names two groups from `GROUP`, each with at least 2 samples.

## Command-line scripts

Each step is an Rscript in `bin/` that you can also run outside Nextflow (with `SPADE_HOME` pointing at a
SPADE checkout). Options follow RNAseqERCC v1.1 (`stdRNAseq.ERCC.Main.R` / `stdRNAseq.ERCC.Reporter.R`):

| Pipeline parameter | Script option | RNAseqERCC option |
|---|---|---|
| `--counts` | `spade_de.R -i/--rawCount` | `-i/--rawCount` |
| `--meta` | `spade_de.R -m/--metaData` | `-m/--metaData` |
| `--comparisons` | `spade_de.R -L/--compList` | `-L/--compList` |
| `--prefix` | `spade_de.R -o/--outPrefix` | `-o/--outPrefix` |
| `--ercc` | `spade_de.R -E/--ERCC` | `-E/--ERCC` |
| `--filter_noncoding` | `spade_de.R -f/--filter` | `-f/--filter` (wrapper `-f` disables) |
| `--cutoff_count` | `spade_de.R -C/--countCutoff` | `-C/--countCutoff` |
| `--cutoff_cpm` | `spade_de.R -c/--cpmCutoff` | `-c/--cpmCutoff` |
| `--cutoff_log2fc`, `--cutoff_pvalue`, `--cutoff_fdr` | report config `CUTOFF_LOG2FC`, `CUTOFF_PVALUE`, `CUTOFF_FDR` | `-F`, `-P`, `-q` / config.yaml |
| `--topn_deg` | report config `NUM_DEG_TOPN` | `-n` / `NUM_DEG_TOPN` |
| `--num_gene_pca` | report config `NUM_GENE_PCA` | `NUM_GENE_PCA` |
| `--volcano` | `spade_report.R -V/--volcano` | wrapper `-V` |
| `--mapping_report` (CAB) | `spade_de.R -r/--report` | `-r/--report` |

Differences from RNAseqERCC v1.1, all from SPADE v3:

- `--ERCC FALSE` really uses TMM normalization (RNAseqERCC still passed ERCC counts to the DE step);
- with no ERCC- rows the run switches to TMM normalization and reports it (RNAseqERCC stopped with an error);
- ERCC normalization factors are used as computed by `edgeR::calcNormFactors()` (RNAseqERCC rescaled them to a
  geometric mean of 1; on SPADE example 3 the DE tables are identical either way);
- `.rnk` files: `_t.rnk` and `_sign.log2FCxNegLog10Pval.rnk` (no `_log2FC.rnk`), with Ensembl IDs when available;
- extra QC outputs (`QC/`: library size, RLE, PCA coordinates, ERCC QC bundle) and GSEA `.cls`/TMM files;
- the results object is `<prefix>_SPADE_Results.rds` (same layout as `*_RNAseqERCC_Results.rds`, plus `settings`).

## Enrichment

**fgsea.** Genes are ranked by the moderated t statistic (`--fgsea_rank_metric t`; or `logPval`,
`log2FC`) using gene symbols, and tested against MSigDB GMT files of release `--msigdb_version` (2026.2):

| `--fgsea_gene_sets default` | Human (`2026.2.Hs`) | Mouse (`2026.2.Mm`) |
|---|---|---|
| GO_MF, GO_BP, GO_CC | `c5.go.mf`, `c5.go.bp`, `c5.go.cc` | `m5.go.mf`, `m5.go.bp`, `m5.go.cc` |
| CellType | `c8.all` | `m8.all` |
| CanonicalPathway | `c2.cp` | `m2.cp` |
| CGP | `c2.cgp` | `m2.cgp` |
| Hallmark | `h.all` | `mh.all` |
| OncogenicSignatures / TumorPhenotype | `c6.all` | `m5.mpo.mpt` |

Other names in [assets/msigdb_gene_sets.csv](../assets/msigdb_gene_sets.csv) (Reactome, WikiPathways, ...) or any
MSigDB file stem (e.g. `c2.cp.reactome`) can be listed: `--fgsea_gene_sets Hallmark,c2.cp.reactome`.
Add your own GMT files (gene symbols) with `--custom_gmt a.gmt,b.gmt`.

The GMT files are downloaded from the Broad Institute when the pipeline starts. If the machine running Nextflow
has no internet, download them once and use `--msigdb_dir`:

```bash
v=2026.2; sp=Hs   # Mm for mouse
for s in c5.go.mf c5.go.bp c5.go.cc c8.all c2.cp c2.cgp h.all c6.all; do
  wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/$v.$sp/$s.v$v.$sp.symbols.gmt
done
```

**Enrichr** (`--run_enrichr true`) sends the top `--enrichr_top_genes` up- and down-regulated genes to the
Enrichr web service, so compute nodes need outbound HTTPS. A failed Enrichr task is ignored (the report is still
made). `--enrichr_dbs standard` picks the newest Hallmark, GO BP/MF/CC, KEGG (human or mouse) and Reactome libraries.

## CAB internal options

For St. Jude CAB analyses of AutoMapper output. They need `--cab_internal`, which also switches the report to
the CAB branding (CAB logo, wiki links, project space paths, acknowledgements):

```bash
# one or more SRMs from the AutoMapper directory (several are merged, as stdRNAseq.ERCC.sh does)
nextflow run hongjianjin/SPADE-nf -profile stjude --cab_internal \
    --automapper_dir /path/to/automapper --srm SRM1_SRM2 \
    --meta sample_meta.txt --comparisons comparisons.lst --species human --outdir results \
    --pi_name "..." --project_name "..." --lead_analysts "..." --mapping_qc good

# or a sheet with one row per SRM: srm,counts,report,tpm,fpkm
nextflow run hongjianjin/SPADE-nf -profile stjude --cab_internal --input srms.csv ...
```

`--mapping_report <SRM>.report` adds the mapping QC tables when `--counts` is used. With several SRMs, count
matrices are inner-joined on geneID (same genes and annotation required; sample IDs must be unique) and the
output prefix is `mergedSRMs_<n>`.

## St. Jude HPC

```bash
module load conda3/202402
conda create -n nextflow -c conda-forge -c bioconda nextflow=24.10.5   # once
source $(conda info --base)/etc/profile.d/conda.sh; conda activate nextflow
export NXF_SINGULARITY_CACHEDIR=/research/groups/cab/projects/Control/common/cab_transcriptomics/dev/singularity_cache
export SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR/.singularity

bsub -q standard -n 1 -M 8000 -W 72:00 -J spade-nf -o spade-nf.%J.log \
    "nextflow run hongjianjin/SPADE-nf -r v3.0.1 -profile stjude --counts ... --outdir results -resume"
```

## Troubleshooting

- **`unable to load shared object '/home/<user>/R/...'`**: a host `R_LIBS` leaked into the container; the pipeline
  resets it, but check `singularity exec --cleanenv docker://ghcr.io/hongjianjin/spade:3.0.1 Rscript -e 'library(limma)'`.
- **Input validation failed**: read `results/input/<prefix>_validation.tsv` and the task's `.command.log`.
- **fgsea downloads fail**: use `--msigdb_dir` (above).
- **Enrichr: cannot reach the Enrichr server**: compute nodes have no internet; drop `--run_enrichr`.
- **Stale container after a rebuild under the same tag**: delete
  `$NXF_SINGULARITY_CACHEDIR/ghcr.io-hongjianjin-spade-<tag>.img`.
