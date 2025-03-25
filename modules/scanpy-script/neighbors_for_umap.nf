process NEIGHBORS_FOR_UMAP {
    container params.scanpy_scripts_container


    input:
        tuple path(anndata), val(n_neighbors)
        val representation
    output:
        path "neighbors_${n_neighbors}.h5ad"
    script:
    """
        export PYTHONIOENCODING='utf-8'
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