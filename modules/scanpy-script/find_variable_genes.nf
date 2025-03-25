process FIND_VARIABLE_GENES {
    container params.scanpy_scripts_container

    input:
        path anndata
        val batch_variable

    output:
        path 'variable_genes.h5ad'

    script:
    """
        batch_variable_tag=""
        if [[ -n "$batch_variable" ]]; then
            batch_variable_tag="--batch-key $batch_variable"
        fi

        export PYTHONIOENCODING='utf-8'
        scanpy-find-variable-genes \
        --flavor 'seurat' \
        --mean-limits 0.0125 1000000000.0 \
        --disp-limits 0.5 50.0 \
        --span 0.3 \
        --n-bins '20' \
        \$batch_variable_tag \
        --input-format 'anndata' \
        $anndata \
        --show-obj stdout \
        --output-format anndata 'variable_genes.h5ad'
    """
}