// Calls the Enrichr web API: needs outbound HTTPS from the compute node.
process ENRICHR {
    tag "${meta.id}"
    label 'process_single'
    label 'error_ignore'

    container "${params.spade_container}"

    input:
    tuple val(meta), path(de_table)
    val species

    output:
    tuple val(meta), path("${prefix}_enrichr.tsv"), emit: results
    path "${prefix}_*.pdf"                         , emit: plots, optional: true
    path "enrichR/*.xlsx"                          , emit: xlsx , optional: true
    path "${prefix}_enrichr_input_*.txt"           , emit: genes, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    spade_enrichr.R \\
        --de ${de_table} \\
        --comparison ${meta.id} \\
        --species ${species} \\
        --outPrefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        enrichr: \$(Rscript -e 'cat(as.character(packageVersion("enrichR")))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'comparison\\tdirection\\tlibrary\\tTerm\\n' > ${prefix}_enrichr.tsv
    printf '"%s":\\n    r-base: stub\\n    enrichr: stub\\n' "${task.process}" > versions.yml
    """
}
