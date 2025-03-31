process HARMONY_BATCH {
    container params.scanpy_scripts_container

    input:
        path anndata
        val batch_variable
    output:
        path 'harmony.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
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
    stub:
    """
        touch harmony.h5ad
    """
}