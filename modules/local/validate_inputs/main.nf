process VALIDATE_INPUTS {
    tag "${meta.id}"
    label 'process_single'

    container "${params.spade_container}"

    input:
    tuple val(meta), path(counts), path(sample_meta), path(comparisons)
    val species

    output:
    tuple val(meta), path("${prefix}_comparisons.lst"), emit: comparisons
    path "${prefix}_validation.tsv"                    , emit: validation
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix        = task.ext.prefix ?: "${meta.id}"
    def comp_arg  = comparisons ? "--compList ${comparisons}" : ''
    """
    spade_validate.R \\
        --rawCount ${counts} \\
        --metaData ${sample_meta} \\
        --species ${species} \\
        ${comp_arg} \\
        --outPrefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        spade: \$(cat \${SPADE_HOME:-/opt/spade}/VERSION 2> /dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ -s "${comparisons}" ]; then cp ${comparisons} ${prefix}_comparisons.lst; else echo "A_vs_B" > ${prefix}_comparisons.lst; fi
    printf 'key\\tvalue\\nstub\\ttrue\\n' > ${prefix}_validation.tsv
    printf '"%s":\\n    r-base: stub\\n    spade: stub\\n' "${task.process}" > versions.yml
    """
}
