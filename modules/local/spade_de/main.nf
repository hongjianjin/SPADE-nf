process SPADE_DE {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.spade_container}"

    input:
    tuple val(meta), path(counts), path(sample_meta), path(comparisons)
    path mapping_report

    output:
    tuple val(meta), path("${prefix}_SPADE_Results.rds"), emit: rds
    path "${prefix}_*_log2CPM.tsv"                       , emit: de_tables
    path "*.rnk"                                          , emit: rnk      , optional: true
    path "${prefix}_all_GSEAinput*"                       , emit: gsea_input, optional: true
    path "${prefix}_ERCC_*"                               , emit: ercc_qc  , optional: true
    path "${prefix}_MappingReport.tsv"                    , emit: mapping  , optional: true
    path "QC"                                             , emit: qc       , optional: true
    path "${prefix}_log.txt"                              , emit: log
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def report   = mapping_report ? "--report ${mapping_report}" : ''
    """
    # The results object is the expensive output; plotting extras after it only log warnings.
    spade_de.R \\
        --rawCount ${counts} \\
        --metaData ${sample_meta} \\
        --compList ${comparisons} \\
        ${report} \\
        --outPrefix ${prefix} \\
        --ncores ${task.cpus} \\
        ${args}

    if [ ! -s ${prefix}_SPADE_Results.rds ]; then
        echo "ERROR: ${prefix}_SPADE_Results.rds was not produced; see ${prefix}_log.txt" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        spade: \$(cat \${SPADE_HOME:-/opt/spade}/VERSION 2> /dev/null || echo unknown)
        limma: \$(Rscript -e 'cat(as.character(packageVersion("limma")))' 2> /dev/null | tail -n 1)
        edger: \$(Rscript -e 'cat(as.character(packageVersion("edgeR")))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_SPADE_Results.rds ${prefix}_log.txt
    while read -r comp; do
        [ -n "\$comp" ] && touch "${prefix}_\${comp}_log2CPM.tsv" "${prefix}_\${comp}_t.rnk"
    done < ${comparisons}
    mkdir -p QC && touch QC/${prefix}_QC_libSize.pdf
    printf '"%s":\\n    r-base: stub\\n    spade: stub\\n' "${task.process}" > versions.yml
    """
}
