process NEIGHBORS {
    container params.scanpy_scripts_container

    input:
        path anndata
        val representation
    output:
        path 'neighbors.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
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
    stub:
    """
        touch neighbors.h5ad
    """
}