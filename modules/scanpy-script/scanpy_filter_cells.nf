process SCANPY_FILTER_CELLS {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val category

    output:
        path 'filtered_cell_anndata.h5ad'

    script:
    """
    n_counts=1500
    if [[ -n "$category" ]]; then
        n_counts=750
    fi

    export PYTHONIOENCODING='utf-8'
	scanpy-filter-cells --gene-name 'gene_symbols' \
        --param 'c:n_counts' \$n_counts 1000000000.0 \
        --param 'c:pct_counts_mito' 0.0 0.35 \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata 'filtered_cell_anndata.h5ad' \
        --export-mtx ./ \
        $category
    """
}