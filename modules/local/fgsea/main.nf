process FGSEA {
    tag "${meta.id}:${gene_set}"
    label 'process_single'

    container "${params.spade_container}"

    input:
    tuple val(meta), path(de_table), val(gene_set), path(gmt)

    output:
    tuple val(meta), path("${prefix}_fgsea.tsv"), emit: results
    path "${prefix}_fgsea_*.pdf"                 , emit: plots, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_${gene_set}"
    """
    spade_fgsea.R \\
        --de ${de_table} \\
        --gmt ${gmt} \\
        --gene_set ${gene_set} \\
        --comparison ${meta.id} \\
        --outPrefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        fgsea: \$(Rscript -e 'cat(as.character(packageVersion("fgsea")))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${gene_set}"
    """
    printf 'comparison\\tgene_set\\tpathway\\tpval\\tpadj\\tNES\\tsize\\n' > ${prefix}_fgsea.tsv
    printf '"%s":\\n    r-base: stub\\n    fgsea: stub\\n' "${task.process}" > versions.yml
    """
}
