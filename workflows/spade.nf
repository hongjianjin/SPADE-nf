/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERGE_SRMS      } from '../modules/local/merge_srms/main'
include { VALIDATE_INPUTS } from '../modules/local/validate_inputs/main'
include { SPADE_DE        } from '../modules/local/spade_de/main'
include { FGSEA           } from '../modules/local/fgsea/main'
include { ENRICHR         } from '../modules/local/enrichr/main'
include { SPADE_REPORT    } from '../modules/local/spade_report/main'
include { runPrefix       } from '../subworkflows/local/utils_spade_pipeline'
include { geneSetList     } from '../subworkflows/local/utils_spade_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPADE {

    take:
    ch_srms     // channel: [ val(meta), path(counts), path(report) | [], path(tpm) | [], path(fpkm) | [] ], one per SRM

    main:
    ch_versions = Channel.empty()
    def prefix    = runPrefix()
    def gene_sets = params.run_fgsea ? geneSetList() : []   // [ [ name, gmt ], ... ]
    def meta0   = [ id: prefix ]
    ch_meta     = Channel.value(file(params.meta, checkIfExists: true))
    ch_comps    = params.comparisons ? file(params.comparisons, checkIfExists: true) : []

    //
    // CAB internal: merge several SRMs into one count matrix (and mapping report)
    //
    ch_srm_list = ch_srms.toList()
    ch_multi = ch_srm_list
        .filter { rows -> rows.size() > 1 }
        .map { rows ->
            [ meta0,
              rows.collect { it[1] },
              rows.collect { it[2] }.findAll { it },
              rows.collect { it[3] }.findAll { it },
              rows.collect { it[4] }.findAll { it } ]
        }
    MERGE_SRMS ( ch_multi )
    ch_versions = ch_versions.mix(MERGE_SRMS.out.versions)

    ch_single = ch_srm_list
        .filter { rows -> rows.size() == 1 }
        .map { rows -> [ meta0, rows[0][1], rows[0][2] ?: [] ] }
    ch_counts_report = ch_single.mix(
        MERGE_SRMS.out.counts
            .join(MERGE_SRMS.out.report, remainder: true)
            .map { meta, counts, report -> [ meta, counts, report ?: [] ] }
    )
    ch_counts = ch_counts_report.map { meta, counts, report -> [ meta, counts ] }
    ch_report = ch_counts_report.map { meta, counts, report -> report ?: [] }
    if (params.mapping_report) {
        ch_report = Channel.value(file(params.mapping_report, checkIfExists: true))
    }

    //
    // MODULE: check metadata, counts, species and comparisons before the expensive step
    //
    VALIDATE_INPUTS (
        ch_counts.combine(ch_meta).map { meta, counts, sm -> [ meta, counts, sm, ch_comps ] },
        params.species
    )
    ch_versions = ch_versions.mix(VALIDATE_INPUTS.out.versions)

    //
    // MODULE: QC, ERCC QC and limma-voom DE for all comparisons (one cohort-level job)
    //
    SPADE_DE (
        ch_counts.combine(ch_meta).join(VALIDATE_INPUTS.out.comparisons)
            .map { meta, counts, sm, comps -> [ meta, counts, sm, comps ] },
        ch_report.ifEmpty([])
    )
    ch_versions = ch_versions.mix(SPADE_DE.out.versions)

    // One DE table per comparison: <prefix>_<comparison>_log2CPM.tsv
    ch_de = SPADE_DE.out.de_tables
        .flatten()
        .map { f ->
            def comp = f.name.replaceFirst("^${prefix}_", '').replaceFirst(/_log2CPM\.tsv$/, '')
            [ [ id: comp ], f ]
        }

    //
    // MODULE: fgsea, one task per comparison x gene set
    //
    ch_enrichment = Channel.empty()
    if (params.run_fgsea && gene_sets) {
        FGSEA (
            ch_de.combine(Channel.fromList(gene_sets))
        )
        ch_enrichment = ch_enrichment.mix(FGSEA.out.results.map { meta, f -> f })
        ch_versions   = ch_versions.mix(FGSEA.out.versions.first())
    }

    //
    // MODULE: Enrichr (web API), one task per comparison
    //
    if (params.run_enrichr) {
        ENRICHR ( ch_de, params.species )
        ch_enrichment = ch_enrichment.mix(ENRICHR.out.results.map { meta, f -> f })
        ch_versions   = ch_versions.mix(ENRICHR.out.versions.first())
    }

    //
    // MODULE: HTML report (cutoffs and project info only change this step)
    //
    ch_config = Channel.of(reportConfig(prefix)).collectFile(name: "${prefix}_report_config.yaml")
    SPADE_REPORT (
        SPADE_DE.out.rds,
        ch_config,
        ch_enrichment.collect().ifEmpty([]),
        file("${projectDir}/assets/report", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SPADE_REPORT.out.versions)

    //
    // Collate software versions
    //
    ch_versions
        .unique()
        .collectFile(
            storeDir : "${params.outdir}/pipeline_info",
            name     : 'spade_software_versions.yml',
            sort     : true,
            newLine  : true
        )

    emit:
    report   = SPADE_REPORT.out.html
    results  = SPADE_DE.out.rds
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// KEY: value file read by spade_report.R (same keys as the RNAseqERCC config)
//
def reportConfig(prefix) {
    def fmt = { v -> v == null ? '' : v.toString() }
    def lines = [
        SRM               : prefix,
        PIPELINE          : "${workflow.manifest.name} v${workflow.manifest.version} (SPADE ${params.spade_version})",
        PROJECT_NAME      : params.project_name,
        PI_NAME           : params.pi_name,
        PROJECT_LEAD_NAME : params.project_lead,
        LEAD_ANALYSTS     : params.lead_analysts,
        GROUP_LEAD        : params.group_lead,
        CONTACT_EMAIL     : params.contact_email,
        DEPARTMENT        : params.department,
        TASK_ID           : params.task_id,
        START_DATE        : params.start_date,
        COMPLETION_DATE   : new java.util.Date().format('yyyy-MM-dd'),
        MAPPING_QC        : params.mapping_qc,
        SPECIES           : params.species,
        OUTPUT_DIR        : params.outdir,
        COUNT_FILE        : params.counts ? file(params.counts).name : 'merged SRM count matrix',
        META_FILE         : file(params.meta).name,
        CUTOFF_COUNT      : params.cutoff_count,
        CUTOFF_CPM        : params.cutoff_cpm,
        NUM_GENE_PCA      : params.num_gene_pca,
        NUM_DEG_TOPN      : params.topn_deg,
        CUTOFF_PVALUE     : params.cutoff_pvalue,
        CUTOFF_FDR        : params.cutoff_fdr,
        CUTOFF_LOG2FC     : params.cutoff_log2fc,
        MSIGDB_VERSION    : params.msigdb_version,
        FGSEA_GENE_SETS   : params.fgsea_gene_sets,
        FGSEA_RANK_METRIC : params.fgsea_rank_metric,
        ENRICHR_DBS       : params.enrichr_dbs
    ]
    // values are single-line; ':' is fine after the first one (ReadYaml splits on the first ':')
    return lines.collect { k, v -> "${k}: ${fmt(v).replaceAll(/[\r\n]+/, ' ')}" }.join('\n') + '\n'
}
