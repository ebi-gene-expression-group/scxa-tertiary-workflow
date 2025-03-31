// a fake process
process SCANPY_MULTIPLET_SCRUBLET1 {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val batch_variable

    output:
        path 'scrublet3.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
        if [ -z "$batch_variable" ]; then
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            $anndata \
            scrublet3.h5ad
        else
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            --batch-key "$batch_variable" \
            $anndata \
            scrublet3.h5ad
        fi
    """
    stub:
    """
        touch scrublet3.h5ad
    """
}

// another fake process
process SCANPY_MULTIPLET_SCRUBLET2 {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val batch_variable

    output:
        path 'scrublet2.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
        if [ -z "$batch_variable" ]; then
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            $anndata \
            scrublet2.h5ad
        else
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            --batch-key "$batch_variable" \
            $anndata \
            scrublet3.h5ad
        fi
    """
    stub:
    """
        touch scrublet2.h5ad
    """
}