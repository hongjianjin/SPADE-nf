// Cheap and tunable: report cutoffs live in the config file, so changing them with -resume
// re-renders the report without re-running the DE step.
process SPADE_REPORT {
    tag "${meta.id}"
    label 'process_low'

    container "${params.spade_container}"

    input:
    tuple val(meta), path(rds)
    path config
    path enrichment, stageAs: 'enrichment/*'
    path report_assets, stageAs: 'report_assets'

    output:
    path "*_DE_Report_*.html"      , emit: html
    path "${prefix}_DEG_summary.tsv", emit: deg_summary
    path "${prefix}_DEG_counts.tsv" , emit: deg_counts
    path "*_fgsea_results.xlsx"    , emit: fgsea_xlsx, optional: true
    path "volcanoPlots"            , emit: volcano   , optional: true
    path "${prefix}_report_log.txt", emit: log
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\$PWD   # pandoc / rmarkdown want a writable HOME inside the container
    spade_report.R \\
        --config ${config} \\
        --input ${rds} \\
        --outDir . \\
        --enrichment_dir enrichment \\
        --assets_dir report_assets \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))' 2> /dev/null | tail -n 1)
        rmarkdown: \$(Rscript -e 'cat(as.character(packageVersion("rmarkdown")))' 2> /dev/null | tail -n 1)
        pandoc: \$(pandoc --version 2> /dev/null | head -n 1 | sed 's/^pandoc //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch SPADE_DE_Report_${prefix}.html ${prefix}_DEG_summary.tsv ${prefix}_DEG_counts.tsv ${prefix}_report_log.txt
    printf '"%s":\\n    r-base: stub\\n' "${task.process}" > versions.yml
    """
}
