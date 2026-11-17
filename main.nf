#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPADE-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ERCC spike-in aware RNA-seq differential expression and enrichment with SPADE
    https://github.com/hongjianjin/SPADE
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { SPADE                   } from './workflows/spade'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_spade_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_spade_pipeline'

workflow {

    main:
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params
    )

    SPADE (
        PIPELINE_INITIALISATION.out.srms
    )

    PIPELINE_COMPLETION (
        params.outdir
    )
}
