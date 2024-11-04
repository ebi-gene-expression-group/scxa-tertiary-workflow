#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.celltype_field = 'NO_CELLTYPE_FIELD'
params.neighbor_values = ['10', '100', '15', '20', '25', '3', '30', '5', '50']
params.perplexity_values = ['1', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50']
params.resolution_values = ['0.1', '0.3', '0.5', '0.7', '1.0', '2.0', '3.0', '4.0', '5.0']
params.slotname = "louvain_resolution"
params.clustering_slotname = params.resolution_values.collect { params.slotname + "_" + it }
params.merged_group_slotname = params.clustering_slotname + params.celltype_field


log.info """
===============================
WORKFLOW PARAMETER VALUES
===============================
celltype_field: ${params.celltype_field}
neighbor_values: ${params.neighbor_values}
perplexity_values: ${params.perplexity_values}
resolution_values: ${params.resolution_values}
slotname: ${params.slotname}
clustering_slotname: ${params.clustering_slotname}
merged_group_slotname: ${params.merged_group_slotname}
===============================
"""

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
        join -a 1 -t \$'\t' -o 0,1.2,2.2 sorted_gene.txt sorted_genemeta.txt | cut -f1,3 > merged_genemeta.tsv
    """
}

process scanpy_read_10x {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'
    
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
        ln -s $genes genes.tsv
        #ln -s $barcodes barcodes.tsv
        
        scanpy-read-10x --input-10x-mtx ./ \
        --var-names 'gene_ids' \
        --extra-obs $cellmeta \
        --extra-var $genemeta \
        --show-obj stdout \
        --output-format anndata \
        'anndata.h5ad'
    """
}

process scanpy_filter_cells {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'
    
    input:
        path anndata
        path genes

    output:
        path 'filtered_cell_anndata.h5ad'

    script:
    """
        scanpy-filter-cells --gene-name 'gene_symbols' \
        --param 'c:n_counts' 750.0 1000000000.0 \
        --param 'c:pct_counts_mito' 0.0 0.35 \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata 'filtered_cell_anndata.h5ad' \
        --export-mtx ./
    """
}

process scanpy_filter_genes {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        path anndata
        path genes

    output:
        path 'filtered_gene_anndata.h5ad'

    script:
    """
        scanpy-filter-genes \
        --param 'g:n_cells' 3.0 1000000000.0 \
        --subset 'g:index' \
        $genes \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata \
        'filtered_gene_anndata.h5ad' \
        --export-mtx ./
    """
}

process normalise_data {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        path anndata

    output:
        path 'normalised_anndata.h5ad'

    script:
    """
        scanpy-normalise-data \
        --no-log-transform \
        --normalize-to '1000000.0' \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata \
        'normalised_anndata.h5ad' \
        --export-mtx ./
    """
}

process normalise_internal_data {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'
    
    input:
        path anndata

    output:
        path 'normalised_internal_anndata.h5ad'

    script:
    """
        scanpy-normalise-data \
        --normalize-to '1000000.0' \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata \
        'normalised_internal_anndata.h5ad' 
    """
}

process find_variable_genes {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        path anndata
        val batch_variable

    output:
        path 'variable_genes.h5ad'

    script:
    """
        batch_variable_tag=""
        if [[ -n "$batch_variable" ]]; then
            batch_variable_tag="--batch-key $batch_variable"
        fi


        scanpy-find-variable-genes \
        --flavor 'seurat' \
        --mean-limits 0.0125 1000000000.0 \
        --disp-limits 0.5 50.0 \
        --span 0.3 \
        --n-bins '20' \
        \$batch_variable_tag \
        --input-format 'anndata' \
        $anndata \
        --show-obj stdout \
        --output-format anndata 'variable_genes.h5ad'
    """
}

process run_pca {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        path anndata

    output:
        path 'PCA.h5ad'

    script:
    """
        scanpy-run-pca \
        --no-zero-center \
        --svd-solver 'arpack' \
        --random-state '1234' \
        --input-format 'anndata' \
        $anndata \
        --show-obj stdout \
        --output-format anndata \
        'PCA.h5ad'
    """
}

process harmony_batch {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        path anndata
        val batch_variable
    output:
        path 'harmony.h5ad'

    script:
    """
        if [[ -n "$batch_variable" ]]; then
            scanpy-integrate harmony \
            --batch-key $batch_variable \
            --basis 'X_pca' \
            --adjusted-basis 'X_pca_harmony' \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            'harmony.h5ad'
        else
            echo "No batch variables passed, simply passing original input as output unchanged."

            cp $anndata 'harmony.h5ad'
        fi

    """
}

process neighbours {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        path anndata
        val pca_param
    output:
        path 'neighbours.h5ad'

    script:
    """
        scanpy-neighbors \
        --n-neighbors 15 \
        --method 'umap' \
        --metric 'euclidean' \
        --random-state '0' \
        --use-rep $pca_param \
        --n-pcs '50' \
        --input-format 'anndata' \
        $anndata \
        --show-obj stdout \
        --output-format anndata \
        'neighbours.h5ad'

    """
}

process neighbours_for_umap {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    memory { 4.GB * task.attempt }
    maxRetries 3

    input:
        tuple path(anndata), val(n_neighbours)
        val pca_param
    output:
        path "neighbours_${n_neighbours}.h5ad"
    script:
    """
        scanpy-neighbors \
            --n-neighbors $n_neighbours \
            --method 'umap' \
            --metric 'euclidean' \
            --random-state '0' \
            --use-rep $pca_param \
            --n-pcs '50' \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            'neighbours_${n_neighbours}.h5ad'

    """
}

process find_clusters {
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'

    input:
        tuple path(anndata), val(resolution)
    output:
        path "clusters_${resolution}.h5ad"
    script:
    """
        scanpy-find-cluster louvain \
        --neighbors-key 'neighbors' \
        --key-added 'louvain_resolution_${resolution}' \
        --resolution ${resolution} \
        --random-state '1234' \
        --directed \
        --export-cluster output.tsv \
        --input-format 'anndata' \
        $anndata \
        --show-obj stdout \
        --output-format anndata \
        'clusters_${resolution}.h5ad'
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
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'
    
    input:
        path anndata
    output:
        path 'umap_*.h5ad'
    script:
    """
            scanpy-run-umap \
            --neighbors-key 'neighbors_$anndata' \
            --key-added 'neighbors_$anndata' \
            --export-embedding embeddings.tsv \
            --n-components 2 \
            --min-dist 0.5 \
            --spread 1.0 \
            --alpha 1.0 \
            --gamma 1.0 \
            --negative-sample-rate 5 \
            --random-state 0 \
            --init-pos 'spectral' \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            'umap_$anndata.h5ad'  
            # Not sure if following is needed
            # && mv 'embeddings_neighbors_n_neighbors_100.tsv' embeddings.tsv

    """
}

process run_tsne {
    //errorStrategy 'ignore'
    
    container 'quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0'
    
    input:
        tuple path(anndata), val(perplexity_values)
        val pca_param
    output:
        path "tsne_${perplexity_values}.h5ad"
    script:
    """
            scanpy-run-tsne \
            --use-rep $pca_param \
            --export-embedding embeddings.tsv \
            --perplexity $perplexity_values \
            --key-added 'perplexity_$perplexity_values' \
            --early-exaggeration '12.0' \
            --learning-rate '400.0' \
            --no-fast-tsne \
            --random-state 1234  \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            'tsne_${perplexity_values}.h5ad'
            # Not sure if following is needed
            # && mv 'embeddings_perplexity_1.tsv' embeddings.tsv
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
    batch_variable = Channel.value('')
    neighbors_ch = channel.fromList(params.neighbor_values)
    perplexity_ch = channel.fromList(params.perplexity_values)
    resolution_ch = channel.fromList(params.resolution_values)

    Column_rearrange_1(
        genemeta, 
        "gene_id"
    )
    Column_rearrange_2(
        genemeta, 
        "gene_id", 
        "gene_name"
    )
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
        scanpy_read_10x.out,
        Column_rearrange_1.out[0]
    )
    normalise_data(
        scanpy_filter_cells.out
    )
    normalise_internal_data(
        scanpy_filter_cells.out
    )
    find_variable_genes(
        normalise_internal_data.out,
        batch_variable
    )
    run_pca(
        find_variable_genes.out
    )
    harmony_batch(
        run_pca.out,
        batch_variable
    )
    neighbours(
        harmony_batch.out,
        pca_param
    )
    neighbours_for_umap(
        harmony_batch.out.combine(neighbors_ch),
        pca_param
    )
    TNSEs_ch = run_tsne(
        harmony_batch.out.combine(perplexity_ch),
        pca_param
    )
    //TNSEs_ch
    //    .filter { it.exitStatus == 0 }

    UMAPs_ch = run_umap(
        neighbours_for_umap.out.flatten()
    )
    //UMAPs_ch
   //     .filter { it.exitStatus == 0 }
    find_clusters(
        neighbours.out.combine(resolution_ch)
    )
}
