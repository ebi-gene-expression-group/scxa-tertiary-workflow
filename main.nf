#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.scanpy_scripts_container = "quay.io/biocontainers/scanpy-scripts:1.1.6--pypyhdfd78af_0"
params.technology = "plate"
params.batch_variable = ""
params.representation = "X_pca"
params.dir_path = "."
params.result_dir_path = params.output_path ?: params.dir_path + "/results"
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
    container params.scanpy_scripts_container
    
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

process scanpy_multiplet_scrublet {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val batch_variable

    output:
        path 'scrublet.h5ad'

    script:
    """
        if [ -z "$batch_variable" ]; then
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            $anndata \
            scrublet.h5ad
        else
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            --batch-key "$batch_variable" \
            $anndata \
            scrublet.h5ad
        fi
    """
}

process scanpy_plot_scrublet {
    publishDir params.result_dir_path, mode: 'copy', pattern: '(scrublet.png)'
    container params.scanpy_scripts_container
    
    input:
        path anndata

    output:
        path 'scrublet.png'

    script:
    """
        scanpy-cli plot scrublet \
        --input-format "anndata" \
        --scale-hist-obs "linear" \
        --scale-hist-sim "linear" \
        $anndata \
        scrublet.png
    """
}

process scanpy_filter_cells {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val category

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
        --export-mtx ./ \
        $category
    """
}

process scanpy_filter_genes {
    publishDir "${params.result_dir_path}/matrices/raw_filtered", mode: 'copy', pattern: 'matrix.mtx'
    publishDir "${params.result_dir_path}/matrices/raw_filtered", mode: 'copy', pattern: 'barcodes.tsv'
    publishDir "${params.result_dir_path}/matrices/raw_filtered", mode: 'copy', pattern: 'genes.tsv'
    container params.scanpy_scripts_container

    input:
        path anndata
        path genes

    output:
        path 'filtered_gene_anndata.h5ad'
        path 'matrix.mtx'
        path 'barcodes.tsv'
        path 'genes.tsv'

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
    publishDir "${params.result_dir_path}/matrices/filtered_normalised", mode: 'copy', pattern: 'matrix.mtx'
    publishDir "${params.result_dir_path}/matrices/filtered_normalised", mode: 'copy', pattern: 'barcodes.tsv'
    publishDir "${params.result_dir_path}/matrices/filtered_normalised", mode: 'copy', pattern: 'genes.tsv'
    container params.scanpy_scripts_container

    input:
        path anndata

    output:
	path 'normalised_anndata.h5ad'
	path 'matrix.mtx'
	path 'barcodes.tsv'
	path 'genes.tsv'

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
    container params.scanpy_scripts_container
    
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
    container params.scanpy_scripts_container

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
    container params.scanpy_scripts_container

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
    container params.scanpy_scripts_container

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

process neighbors {
    container params.scanpy_scripts_container

    input:
        path anndata
        val representation
    output:
        path 'neighbors.h5ad'

    script:
    """
        scanpy-neighbors \
        --n-neighbors 15 \
        --method 'umap' \
        --metric 'euclidean' \
        --random-state '0' \
        --use-rep $representation \
        --n-pcs '50' \
        --input-format 'anndata' \
        $anndata \
        --show-obj stdout \
        --output-format anndata \
        'neighbors.h5ad'

    """
}

process neighbors_for_umap {
    container params.scanpy_scripts_container


    input:
        tuple path(anndata), val(n_neighbors)
        val representation
    output:
        path "neighbors_${n_neighbors}.h5ad"
    script:
    """
        scanpy-neighbors \
            --n-neighbors $n_neighbors \
            --key-added 'neighbors_n_neighbors_${n_neighbors}' \
            --method 'umap' \
            --metric 'euclidean' \
            --random-state '0' \
            --use-rep $representation \
            --n-pcs '50' \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            'neighbors_${n_neighbors}.h5ad'

    """
}

process find_clusters {
    publishDir "${params.result_dir_path}/clusters", mode: 'copy', pattern: 'clusters_*.tsv'
    container params.scanpy_scripts_container

    input:
        tuple path(anndata), val(resolution)
    output:
        path "clusters_${resolution}.h5ad"
	path "clusters_${resolution}.tsv"

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
	
	mv 'output.tsv' 'clusters_${resolution}.tsv'
    """
}

process restore_unscaled {
    container params.scanpy_scripts_container

    input:
	tuple path(anndata), path(normalise_internal_data)

    output:
	path "restore_unscaled_output_${anndata}.h5"

    script:
    """
	ln -s $anndata input.h5
	ln -s $normalise_internal_data r_source.h5
	python ${projectDir}/scripts/restore_unscaled.py
	mv output.h5 'restore_unscaled_output_${anndata}.h5'
    """
}

process find_markers {
    publishDir "${params.result_dir_path}/markers", mode: 'copy', pattern: 'markers_*.tsv'
    errorStrategy 'ignore'
    container params.scanpy_scripts_container

    input:
	tuple path(anndata), val(merged_group_slotname)

    output:
	path "markers_${merged_group_slotname}.h5ad"
	path "markers_${merged_group_slotname}.tsv"

    script:
    """
	scanpy-find-markers \
	--save 'markers_${merged_group_slotname}.tsv' \
	--n-genes '100' \
	--groupby '${merged_group_slotname}' \
	--key-added 'markers_${merged_group_slotname}' \
	--method 'wilcoxon' \
	--use-raw  \
	--reference 'rest' \
	--filter-params 'min_in_group_fraction:0.0,max_out_group_fraction:1.0,min_fold_change:1.0'  \
	--input-format 'anndata' \
	$anndata  \
	--show-obj stdout \
	--output-format anndata \
	'markers_${merged_group_slotname}.h5ad'
    """
}

process run_umap {
    publishDir "${params.result_dir_path}/umap", mode: 'copy', pattern: 'umap_n_neighbors_*.tsv'

    errorStrategy 'ignore'

    container params.scanpy_scripts_container
    
    input:
        path anndata

    output:
        path "umap_*.h5ad"
	path "embeddings_neighbors_neighbors_*.tsv"

    script:
    """
	VAR="$anndata"
	n_number="\${VAR%.h5ad}"
	echo \$n_number
	scanpy-run-umap \
            --neighbors-key "neighbors_n_\${n_number}" \
            --key-added "neighbors_\${n_number}" \
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
            "umap_\${n_number}.h5ad" \
            && mv "embeddings_neighbors_neighbors_\${n_number}.tsv" umap_n_neighbors_\${n_number}.tsv

    """
}

process run_tsne {
    publishDir "${params.result_dir_path}/tsne", mode: 'copy', pattern: 'tsne_perplexity_*\\.tsv'

    errorStrategy 'ignore'
    
    container params.scanpy_scripts_container
    
    input:
        tuple path(anndata), val(perplexity_values)
        val representation

    output:
        path "tsne_${perplexity_values}.h5ad"
	path "embeddings_perplexity_${perplexity_values}.tsv"

    script:
    """
            scanpy-run-tsne \
            --use-rep $representation \
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
            'tsne_${perplexity_values}.h5ad' \
            && mv 'embeddings_perplexity_${perplexity_values}.tsv' 'tsne_perplexity_${perplexity_values}.tsv'
    """
}

process make_project_file {
    publishDir params.result_dir_path, mode: 'copy'

    container params.scanpy_scripts_container

    input:
        path neighbors
        path scanpy_read_10x
        path filter_genes
        path normalise_data
        path find_markers
        path TNSEs_mix_UMAPs
    output:
        path "project.h5ad"
    script:
    """
        ln -s $neighbors input.h5
        ln -s $scanpy_read_10x r_source.h5
        ln -s '$filter_genes' x_source_0.h5
        ln -s '$normalise_data' x_source_1.h5
        count=0
        for i in $find_markers
        do
                ln -sf "\${i}" obs_source_\${count}.h5
                ln -sf "\${i}" uns_source_\${count}.h5
                count=\$((count + 1))
                echo "\${count}"
        done
        count=0
        for i in $TNSEs_mix_UMAPs
        do
                ln -sf "\${i}" embedding_source_\${count}.h5
                count=\$((count + 1))
                echo "\${count}"
        done
        python ${projectDir}/scripts/final_project.py
	mv output.h5 project.h5ad
    """
}

workflow {

    // Create input channel (single file via CLI parameter)
    genemeta = Channel.fromPath("${params.dir_path}/genes_metadata.tsv")
    genes = Channel.fromPath("${params.dir_path}/genes.tsv")
    barcodes = Channel.fromPath("${params.dir_path}/barcodes.tsv")
    matrix = Channel.fromPath("${params.dir_path}/matrix.mtx")
    cellmeta = Channel.fromPath("${params.dir_path}/cell_metadata.tsv")
    neighbors_ch = channel.fromList(params.neighbor_values)
    perplexity_ch = channel.fromList(params.perplexity_values)
    resolution_ch = channel.fromList(params.resolution_values)
    merged_group_slotname_ch = Channel.fromList(params.merged_group_slotname)

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

    if ( params.technology == "droplet" ) {
        SCRUBLET_ch = scanpy_multiplet_scrublet(
            scanpy_read_10x.out,
            params.batch_variable
        )
        scanpy_plot_scrublet(
            SCRUBLET_ch
        )
        scanpy_filter_cells(
            SCRUBLET_ch,
            "--category predicted_doublet False"
        )
    }
    else {
        scanpy_filter_cells(
            scanpy_read_10x.out,
            ""
        )
    }

    scanpy_filter_genes(
        scanpy_filter_cells.out,
        Column_rearrange_1.out[0]
    )
    normalise_data(
        scanpy_filter_genes.out[0]
    )
    normalise_internal_data(
        scanpy_filter_genes.out[0]
    )
    find_variable_genes(
        normalise_internal_data.out,
        params.batch_variable
    )
    run_pca(
        find_variable_genes.out
    )
    harmony_batch(
        run_pca.out,
        params.batch_variable
    )
    neighbors(
        harmony_batch.out,
        params.representation
    )
    neighbors_for_umap(
        harmony_batch.out.combine(neighbors_ch),
        params.representation
    )
    TNSEs_ch = run_tsne(
        harmony_batch.out.combine(perplexity_ch),
        params.representation
    )[0]
    //TNSEs_ch
    //    .filter { it.exitStatus == 0 }

    UMAPs_ch = run_umap(
        neighbors_for_umap.out.flatten()
    )[0]
    //UMAPs_ch
   //     .filter { it.exitStatus == 0 }
    find_clusters(
        neighbors.out.combine(resolution_ch)
    )

    // Combine the outputs of find_clusters and neighbors processes
    combined_outputs = find_clusters.out[0].mix(neighbors.out)

    if ( params.technology == "droplet" ) {
        restore_unscaled (
	    combined_outputs.combine(normalise_internal_data.out)
    	)
	restore_unscaled_files = restore_unscaled.out.map { file ->
	    // Extract the sample number from the file name
	    def sampleNumber = file.baseName.replaceFirst('restore_unscaled_output_', '').replaceFirst('clusters', params.slotname).replaceFirst('neighbors',params.celltype_field).replaceFirst('.h5ad','')
	    [file, sampleNumber] // Create a tuple with sample number and file
	}
	find_markers(
	    restore_unscaled_files
    	)
    }
    else {
	processed_files = combined_outputs.map { file ->
         // Extract the sample number from the file name
         def sampleNumber = file.baseName.replaceFirst('clusters', params.slotname).replaceFirst('neighbors',params.celltype_field)
         [file, sampleNumber] // Create a tuple with sample number and file
	}
        find_markers(
	    processed_files
    	)
    }
    make_project_file(
	neighbors.out,
	scanpy_read_10x.out,
	scanpy_filter_genes.out[0],
	normalise_data.out[0],
	find_markers.out[0].collect(),
	TNSEs_ch.mix(UMAPs_ch).collect()
    )
}
