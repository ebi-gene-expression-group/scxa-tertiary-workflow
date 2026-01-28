#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """
===============================
WORKFLOW PARAMETER VALUES
===============================
EXP dir path: ${params.dir_path}
Selected technology: ${params.technology}
Results results_dir_path: ${params.result_dir_path}
celltype_field: ${params.celltype_field}
neighbor_values: ${params.neighbor_values}
perplexity_values: ${params.perplexity_values}
resolution_values: ${params.resolution_values}
slotname: ${params.slotname}
clustering_slotname: ${params.clustering_slotname}
merged_group_slotname: ${params.merged_group_slotname}
batch_field: ${params.batch_field}
representation: ${params.representation}
===============================
"""

include { SCXA_TERTIARY } from "${projectDir}/workflows/scxa_tertiary.nf"

workflow {
    SCXA_TERTIARY()
}
