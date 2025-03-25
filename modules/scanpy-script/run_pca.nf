process RUN_PCA {
    container params.scanpy_scripts_container

    input:
        path anndata

    output:
        path 'PCA.h5ad'

    script:
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
}