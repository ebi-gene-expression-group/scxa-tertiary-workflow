process HARMONY_BATCH {
    container params.scanpy_scripts_container

    input:
        path anndata
        val batch_field
    output:
        path 'harmony.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        ANNDATA=${WorkflowParamValidator.shellQuote(anndata)}
        BATCH_FIELD=${WorkflowParamValidator.shellQuote(batch_field)}
        export PYTHONIOENCODING='utf-8'
        if [[ -n "\$BATCH_FIELD" ]]; then
            scanpy-integrate harmony \
            --batch-key "\$BATCH_FIELD" \
            --basis 'X_pca' \
            --adjusted-basis 'X_pca_harmony' \
            --input-format 'anndata' \
            "\$ANNDATA" \
            --show-obj stdout \
            --output-format anndata \
            'harmony.h5ad'
        else
            echo "No batch variables passed, simply passing original input as output unchanged."

            cp "\$ANNDATA" 'harmony.h5ad'
        fi

    """
    stub:
    """
        touch harmony.h5ad
    """
}
