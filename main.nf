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
      if [[ -z "\$col1_num" || -z "\$col2_num" ]]; then
          echo "Error: Column '$col1' or '$col2' not found in $genemeta" >&2
          exit 1
      fi
  
      # Extract the gene_id column (without the header)
      tail -n +2 "$genemeta" | cut -f\$col1_num,\$col2_num > filtered_genemeta.txt
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
      path params.output

    script:
    """
        # Sort both files by the first column for join compatibility
        sort -k1,1 "$gene" > sorted_gene.txt
        sort -k1,1 filtered_genemeta.txt > sorted_genemeta.txt
        
        # Perform a left join to keep all data from gene file
        join -a 1 -e 'NA' -t '\t' sorted_gene.txt sorted_genemeta.txt | cut -f1,4 > ${params.output} 
    """
}

process scanpy_read_10x {
    input:

    output:
        path anndata

    conda:
        

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

workflow {

    // Create input channel (single file via CLI parameter)
    genemeta = Channel.fromPath('genemeta_data.txt')
    genes = Channel.fromPath('genes_data.txt')
    barcodes = Channel.fromPath('barcodes_data.txt')
    matrix = Channel.fromPath('matrix_data.txt')
    cellmeta = Channel.fromPath('cellmeta_data.txt')
    pca_param = Channel.value('X_pca')
    celltype_field_param = Channel.value('NO_CELLTYPE_FIELD')
    batch_variable = Channel.value('')
    perplexity_values = Channel.value(['1', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'])
    resolution_values = Channel.value(['0.1', '0.3', '0.5', '0.7', '1.0', '2.0', '3.0', '4.0', '5.0'])


    // Create index file for input BAM file
    Column_rearrange_1(genemeta, "gene_id")
    Column_rearrange_2(genemeta, "gene_id", "gene_name")
}
