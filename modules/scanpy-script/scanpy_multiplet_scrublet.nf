process SCANPY_MULTIPLET_SCRUBLET {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val batch_variable

    output:
        path 'scrublet.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
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
    stub:
    """
        touch scrublet.h5ad
    """
}