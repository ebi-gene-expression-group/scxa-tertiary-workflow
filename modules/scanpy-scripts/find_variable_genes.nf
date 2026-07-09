process FIND_VARIABLE_GENES {
    container params.scanpy_scripts_container

    input:
        path anndata
        val batch_field

    output:
        path 'variable_genes.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        ANNDATA=${WorkflowParamValidator.shellQuote(anndata)}
        BATCH_FIELD=${WorkflowParamValidator.shellQuote(batch_field)}
        batch_field_tag=()
        if [[ -n "\$BATCH_FIELD" ]]; then
            batch_field_tag=(--batch-key "\$BATCH_FIELD")
        fi

        export PYTHONIOENCODING='utf-8'
        scanpy-find-variable-genes \
        --flavor 'seurat' \
        --mean-limits 0.0125 1000000000.0 \
        --disp-limits 0.5 50.0 \
        --span 0.3 \
        --n-bins '20' \
        "\${batch_field_tag[@]}" \
        --input-format 'anndata' \
        "\$ANNDATA" \
        --show-obj stdout \
        --output-format anndata 'variable_genes.h5ad'
    """
    stub:
    """
        touch variable_genes.h5ad
    """
}
