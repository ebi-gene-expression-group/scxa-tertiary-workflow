#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.scanpy_scripts_container = "quay.io/biocontainers/scanpy-scripts:1.1.2--pypyhdfd78af_1"
params.technology = "plate"
params.batch_variable = ""
params.representation = "X_pca"
params.dir_path = "."
params.result_dir_path = params.output_path ?: params.dir_path + "/results"
params.celltype_field = 'NO_CELLTYPE_FIELD'
params.neighbor_values = ['3', '5', '10', '15', '20', '25', '30', '50', '100']
params.perplexity_values = ['1', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50']
params.resolution_values = ['0.1', '0.3', '0.5', '0.7', '1.0', '2.0', '3.0', '4.0', '5.0']
params.slotname = "leiden_resolution"
params.clustering_slotname = params.resolution_values.collect { params.slotname + "_" + it }
params.merged_group_slotname = params.clustering_slotname + [params.celltype_field]

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
batch_variable: ${params.batch_variable}
representation: ${params.representation}
===============================
"""

include { SCXA_TERTIARY } from "${projectDir}/workflows/scxa_tertiary.nf"

workflow {
    SCXA_TERTIARY()
}
