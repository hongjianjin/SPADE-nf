// CAB internal: merge AutoMapper deliverables of several SRMs into one set of files.
process MERGE_SRMS {
    tag "${meta.id}"
    label 'process_single'

    container "${params.spade_container}"

    input:
    tuple val(meta), path(counts, stageAs: 'counts/?/*'), path(reports, stageAs: 'reports/?/*'), path(tpm, stageAs: 'tpm/?/*'), path(fpkm, stageAs: 'fpkm/?/*')

    output:
    tuple val(meta), path("${prefix}_RSEM_gene_count.txt"), emit: counts
    tuple val(meta), path("${prefix}.report")             , emit: report, optional: true
    path "${prefix}_RSEM_gene_TPM.txt"                     , emit: tpm   , optional: true
    path "${prefix}_RSEM_gene_FPKM.txt"                    , emit: fpkm  , optional: true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def joined = { x -> [ x ].flatten().findAll { it }.join(',') }
    def rep  = joined(reports) ? "--reports ${joined(reports)}" : ''
    def t    = joined(tpm)     ? "--tpm ${joined(tpm)}"         : ''
    def f    = joined(fpkm)    ? "--fpkm ${joined(fpkm)}"       : ''
    """
    spade_merge_srms.R \\
        --counts ${joined(counts)} \\
        ${rep} ${t} ${f} \\
        --outPrefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_RSEM_gene_count.txt
    printf '"%s":\\n    r-base: stub\\n' "${task.process}" > versions.yml
    """
}
