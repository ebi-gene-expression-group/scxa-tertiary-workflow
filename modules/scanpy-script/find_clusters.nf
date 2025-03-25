process FIND_CLUSTERS {
    publishDir "${params.result_dir_path}/clusters", mode: 'copy', pattern: 'clusters_resolution_*.tsv'
    container params.scanpy_scripts_container

    input:
        tuple path(anndata), val(resolution)
    output:
        path "clusters_${resolution}.h5ad"
	path "clusters_resolution_${resolution}.tsv"

    script:
    """
        export PYTHONIOENCODING='utf-8'
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
	
	mv 'output.tsv' 'clusters_resolution_${resolution}.tsv'
    """
}