process NORMALISE_INTERNAL_DATA {
    container params.scanpy_scripts_container
    
    input:
        path anndata

    output:
        path 'normalised_internal_anndata.h5ad'

    script:
    """
        export PYTHONIOENCODING='utf-8'
        scanpy-normalise-data \
        --normalize-to '1000000.0' \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata \
        'normalised_internal_anndata.h5ad' 
    """
}