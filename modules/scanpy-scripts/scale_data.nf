process SCALE_DATA {
    container params.scanpy_scripts_container

    input:
        path anndata

    output:
        path 'scaled_anndata.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
        scanpy-scale-data \
        --input-format "anndata" \
        --output-format "anndata" \
        $anndata \
        'scaled_anndata.h5ad'
    """
    stub:
    """
        touch scaled_anndata.h5ad
    """
}