#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Column_rearrange_1: Only keeps the specified columns and removes header
 */
process Column_rearrange_1 {
    // Set the output file
    input:
      path genemeta
      val col

    output:
      path 'filtered_genemeta.txt'

    script:
    """
      # Find the column number of the specified gene_id column name
      col_num=\$(head -n1 "$genemeta" | tr '\\t' '\\n' | grep -n "^$col\$" | cut -d: -f1)
  
      # If column is found, extract it; otherwise, raise an error
      if [[ -z "\$col_num" ]]; then
          echo "Error: Column '$col' not found in $genemeta" >&2
          exit 1
      fi
  
      # Extract the gene_id column (without the header)
      tail -n +2 "$genemeta" | cut -f\$col_num > filtered_genemeta.txt
    """
}

/*
 * Column_rearrange_2: Only keeps the specified columns and removes header
 */
process Column_rearrange_2 {
    // Set the output file
    input:
      path genemeta
      val col1
      val col2

    output:
      path 'filtered_genemeta_2.txt'

    script:
    """
      # Find the column number of the specified gene_id column name
      col_num_1=\$(head -n1 "$genemeta" | tr '\\t' '\\n' | grep -n "^$col1\$" | cut -d: -f1)
      col_num_2=\$(head -n1 "$genemeta" | tr '\\t' '\\n' | grep -n "^$col2\$" | cut -d: -f1)
  
      # If either column is not found, raise an error
      if [[ -z "\$col_num_1" || -z "\$col_num_2" ]]; then
          echo "Error: Column '$col1' or '$col2' not found in $genemeta" >&2
          exit 1
      fi
  
      # Extract the gene_id column (without the header)
      tail -n +2 "$genemeta" | cut -f\$col_num_1,\$col_num_2 > filtered_genemeta_2.txt
    """
}

/*
 * mergeGeneFiles: Merges gene file with genemeta on column 1, and keeps column1 and 4
 */
process mergeGeneFiles {
    input:
      path gene
      path filtered_genemeta

    output:
      path 'merged_genemeta.tsv'

    script:
    """
        # Sort both files by the first column for join compatibility
        sort -k1,1 "$gene" > sorted_gene.txt
        sort -k1,1 "$filtered_genemeta" > sorted_genemeta.txt
        
        # Perform a left join to keep all data from gene file
        join -a 1 -e 'NA' -t '\t' sorted_gene.txt sorted_genemeta.txt | cut -f1,4 > merged_genemeta.tsv
    """
}

process scanpy_read_10x {
    input:
        path matrix
        path genes
        path barcodes
        path cellmeta
        path genemeta

    output:
        path 'anndata.h5ad'

    script:
    """
        #ln -s $matrix matrix.mtx
        #ln -s $genes genes.tsv
        #ln -s $barcodes barcodes.tsv
        
        scanpy-read-10x --input-10x-mtx ./ \
        --var-names 'gene_ids' \
        --extra-obs $cellmeta \
        --extra-var $genemeta  \
        --show-obj stdout \
        --output-format anndata \
        'anndata.h5ad'
    """
}

process scanpy_filter_cells {
    input:
        path anndata

    output:
        path 'filtered_cell_anndata.h5ad'

    script:
    """
        scanpy-filter-cells --gene-name 'gene_symbols' \
        --param 'c:n_counts' 750.0 1000000000.0 \
        --param 'c:pct_counts_mito' 0.0 0.35 \
        --category 'c:predicted_doublet' 'False' \
        --input-format 'anndata' input.h5  \
        --show-obj stdout \
        --output-format anndata 'filtered_cell_anndata.h5ad'  \
        --export-mtx ./
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

workflow {

    // Create input channel (single file via CLI parameter)
    genemeta = Channel.fromPath('genes_metadata.tsv')
    genes = Channel.fromPath('genes.tsv')
    barcodes = Channel.fromPath('barcodes.tsv')
    matrix = Channel.fromPath('matrix.mtx')
    cellmeta = Channel.fromPath('cell_metadata.tsv')
    pca_param = Channel.value('X_pca')
    celltype_field_param = Channel.value('NO_CELLTYPE_FIELD')
    batch_variable = Channel.value('')
    perplexity_values = Channel.value(['1', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'])
    resolution_values = Channel.value(['0.1', '0.3', '0.5', '0.7', '1.0', '2.0', '3.0', '4.0', '5.0'])


    // Create index file for input BAM file
    Column_rearrange_1(genemeta, "gene_id")
    Column_rearrange_2(genemeta, "gene_id", "gene_name")
    mergeGeneFiles(
        genes,
        Column_rearrange_2.out
    )
    scanpy_read_10x(
        matrix,
        mergeGeneFiles.out,
        barcodes,
        cellmeta,
        genemeta
    )  
    scanpy_filter_cells(
        scanpy_read_10x.out
    )  
}
