#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define inputs as channels
Channel.fromPath('genemeta_data.txt').set { genemeta }
Channel.fromPath('genes_data.txt').set { genes }
Channel.fromPath('barcodes_data.txt').set { barcodes }
Channel.fromPath('matrix_data.txt').set { matrix }
Channel.fromPath('cellmeta_data.txt').set { cellmeta }
Channel.value('X_pca').set { pca_param }
Channel.value('NO_CELLTYPE_FIELD').set { celltype_field_param }
Channel.value('').set { batch_variable }
Channel.value(['1', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50']).set { perplexity_values }
Channel.value(['0.1', '0.3', '0.5', '0.7', '1.0', '2.0', '3.0', '4.0', '5.0']).set { resolution_values }

/*
 * Column_rearrange_1: Only keeps the specified columns and removes header
 */
process Column_rearrange_1 {
    input:

    output:

    script:
    """
    """
}

/*
 * Column_rearrange_2: Only keeps the specified columns and removes header
 */
process Column_rearrange_2 {
    // Set the output file
    input:

    output:

    script:
    """
    """
}

/*
 * mergeGeneFiles: Merges gene file with genemeta on column 1, and keeps column1 and 4
 */
process mergeGeneFiles {
    input:

    output:

    script:
    """
    """
}

process scanpy_read_10x {
    input:

    output:

    script:
    """
    """
}

process scanpy_filter_cells {
    input:

    output:

    script:
    """
    """
}

process scanpy_filter_genes {
    input:

    output:

    script:
    """
    """
}

process normalise_data {
    input:

    output:

    script:
    """
    """
}

process normalise_data_internal {
    input:

    output:

    script:
    """
    """
}

process find_variable_genes {
    input:

    output:

    script:
    """
    """
}

process run_pca {
    input:

    output:

    script:
    """
    """
}

process harmony_batch {
    input:

    output:

    script:
    """
    """
}

process neighbours {
    input:

    output:

    script:
    """
    """
}

process neighbours_for_umap {
    input:

    output:

    script:
    """
    """
}

process normalise_data {
    input:

    output:

    script:
    """
    """
}

process find_clusters {
    input:

    output:

    script:
    """
    """
}

process meta_vars {
    input:

    output:

    script:
    """
    """
}

process clustering_slotnames {
    input:

    output:

    script:
    """
    """
}

process merge_group_slotnames {
    input:

    output:

    script:
    """
    """
}

process merge_collections {
    input:

    output:

    script:
    """
    """
}

process build_list {
    input:

    output:

    script:
    """
    """
}

process find_markers {
    input:

    output:

    script:
    """
    """
}

process filtered_cellgroup_markers {
    input:

    output:

    script:
    """
    """
}

process run_umap {
    input:

    output:

    script:
    """
    """
}

process run_tsne {
    input:

    output:

    script:
    """
    """
}

process filter_failed_umap {
    input:

    output:

    script:
    """
    """
}

process filer_failed_tsne {
    input:

    output:

    script:
    """
    """
}

process merge_embeddings {
    input:

    output:

    script:
    """
    """
}



process make_project_file {
    input:

    output:

    script:
    """
    """
}
