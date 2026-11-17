//
// Subworkflows with functionality specific to the SPADE-nf pipeline
//

include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog   } from 'plugin/nf-schema'
include { samplesheetToList  } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version         // boolean: Display version and exit
    validate_params // boolean: Validate parameters against the schema at runtime

    main:

    if (version) {
        log.info "${workflow.manifest.name} ${getWorkflowVersion()}"
        System.exit(0)
    }

    if (validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    checkInputCombination()

    log.info "Count matrices: ${srmRows().collect { it[0].id }.join(', ')}; output prefix: ${runPrefix()}"

    emit:
    srms = Channel.fromList(srmRows()) // channel: [ val(meta), path(counts), report, tpm, fpkm ], one per SRM
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir // path: Path to output directory where results will be published

    main:
    // Keep `workflow` references out of the handler closure: Nextflow runs it with the
    // workflow metadata as delegate, so `workflow.success` there does not resolve.
    workflow.onComplete {
        completionSummary(outdir)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to the troubleshooting section in docs/usage.md"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Count matrix input: one row per SRM ([meta, counts, report, tpm, fpkm])
//
def srmRows() {
    if (params.counts) {
        return [ [ [ id: 'counts' ], file(params.counts, checkIfExists: true), [], [], [] ] ]
    }
    if (params.input) {
        return samplesheetToList(params.input, "${projectDir}/assets/schema_input.json").collect { row ->
            def (meta, counts, report, tpm, fpkm) = row
            [ meta, counts, report ?: [], tpm ?: [], fpkm ?: [] ]
        }
    }
    return automapperRows(params.automapper_dir, params.srm)
}

//
// Output prefix: --prefix, else SPADE (--counts), the SRM name, or mergedSRMs_<n> as in RNAseqERCC
//
def runPrefix() {
    if (params.prefix) return params.prefix.toString()
    if (params.counts) return 'SPADE'
    def n = params.input ? srmRows().size() : params.srm.toString().split('_').findAll { it }.size()
    return n == 1 ? (params.input ? srmRows()[0][0].id : params.srm.toString()) : "mergedSRMs_${n}"
}

//
// Input rules: --counts for everyone; SRM sheet, AutoMapper directory and mapping report
// are CAB internal and need --cab_internal.
//
def checkInputCombination() {
    def cab_inputs = [ 'input', 'automapper_dir', 'srm', 'mapping_report' ].findAll { params[it] }
    if (cab_inputs && !params.cab_internal) {
        error("--${cab_inputs.join(', --')} ${cab_inputs.size() > 1 ? 'are' : 'is a'} CAB internal option${cab_inputs.size() > 1 ? 's' : ''}; add --cab_internal or use --counts.")
    }
    def sources = [ params.counts ? 'counts' : null, params.input ? 'input' : null, params.automapper_dir ? 'automapper_dir' : null ].findAll { it }
    if (sources.size() != 1) {
        error("Give exactly one count source: --counts <matrix>${params.cab_internal ? ', --input <srm sheet> or --automapper_dir <dir> --srm <SRM1_SRM2>' : ''} (got: ${sources ? sources.join(', ') : 'none'}).")
    }
    if (params.automapper_dir && !params.srm) {
        error("--automapper_dir needs --srm (one SRM, or several joined by '_').")
    }
    if (params.species == 'other' && params.run_fgsea && !params.custom_gmt) {
        error("--species other needs --custom_gmt for fgsea (or --run_fgsea false).")
    }
}

//
// CAB internal: find AutoMapper deliverables like stdRNAseq.ERCC.sh does
//   <dir>/<srm>/<srm>_RSEM_count.txt, else the newest <srm>_RSEM_gene_count*.txt (same for TPM / FPKM)
//   <dir>/<srm>/<srm>*.report (newest)
//
def automapperRows(dir, srm) {
    def srms = srm.toString().split('_').findAll { it }
    return srms.collect { s ->
        def d = "${dir}/${s}"
        def pick = { String plain, String glob ->
            def f = file("${d}/${plain}")
            if (f.exists()) return f
            def hits = files("${d}/${glob}").sort { it.name }
            return hits ? hits[-1] : []
        }
        def counts = pick("${s}_RSEM_count.txt", "${s}_RSEM_gene_count*.txt")
        if (!counts) error("No count matrix for SRM ${s} in ${d} (${s}_RSEM_count.txt or ${s}_RSEM_gene_count*.txt)")
        def reports = files("${d}/${s}*.report").sort { it.name }
        [ [ id: s ], counts, reports ? reports[-1] : [],
          pick("${s}_RSEM_TPM.txt", "${s}_RSEM_gene_TPM*.txt"),
          pick("${s}_RSEM_FPKM.txt", "${s}_RSEM_gene_FPKM*.txt") ]
    }
}

//
// fgsea gene sets: [ [ name, gmt ], ... ]
//   default / short names from assets/msigdb_gene_sets.csv, or any MSigDB file stem (e.g. c2.cp.reactome),
//   from --msigdb_dir when given, else downloaded from --msigdb_base_url; plus --custom_gmt files.
//
def geneSetList() {
    def sets = []
    if (params.species in [ 'human', 'mouse' ] && params.fgsea_gene_sets) {
        def sp     = params.species == 'human' ? 'Hs' : 'Mm'
        def col    = params.species
        def table  = file("${projectDir}/assets/msigdb_gene_sets.csv").readLines().findAll { it.trim() && !it.startsWith('#') }
        def header = table[0].split(',', -1) as List
        def rows   = table.drop(1).collect { l -> [ header, l.split(',', -1) as List ].transpose().collectEntries() }
        def wanted = params.fgsea_gene_sets.toString().split(',').collect { it.trim() }.findAll { it }
        if (wanted == [ 'default' ]) {
            wanted = rows.findAll { it['default'] == 'yes' && it[col] }.collect { it.name }
        }
        wanted.each { w ->
            def row  = rows.find { it.name.equalsIgnoreCase(w) }
            def stem = row ? row[col] : w
            if (!stem) {
                log.warn "Gene set ${w} is not available for ${params.species}; skipped."
                return
            }
            def fname = "${stem}.v${params.msigdb_version}.${sp}.symbols.gmt"
            def gmt   = params.msigdb_dir ? file("${params.msigdb_dir}/${fname}", checkIfExists: true)
                                          : file("${params.msigdb_base_url}/${params.msigdb_version}.${sp}/${fname}")
            sets << [ row ? row.name : stem, gmt ]
        }
    }
    if (params.custom_gmt) {
        params.custom_gmt.toString().split(',').collect { it.trim() }.findAll { it }.each { g ->
            def f = file(g, checkIfExists: true)
            sets << [ f.baseName.replaceAll(/[^A-Za-z0-9._-]/, '_'), f ]
        }
    }
    return sets
}

//
// Print a one-line completion summary
//
def completionSummary(outdir) {
    if (workflow.success) {
        log.info "-[${workflow.manifest.name}] Pipeline completed successfully. Results: ${outdir}"
    } else {
        log.error "-[${workflow.manifest.name}] Pipeline completed with errors"
    }
}

//
// Version string for the pipeline
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }
    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }
    return version_string
}
