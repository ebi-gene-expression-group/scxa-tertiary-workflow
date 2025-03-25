process RUN_PCA {
    container params.scanpy_scripts_container

    input:
        path anndata

    output:
        path 'PCA.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
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
    stub:
    """
        touch PCA.h5ad
    """
}