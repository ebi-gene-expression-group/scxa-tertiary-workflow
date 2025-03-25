process NORMALISE_DATA {
    publishDir "${params.result_dir_path}/matrices/filtered_normalised", mode: 'copy', pattern: 'matrix.mtx'
    publishDir "${params.result_dir_path}/matrices/filtered_normalised", mode: 'copy', pattern: 'barcodes.tsv'
    publishDir "${params.result_dir_path}/matrices/filtered_normalised", mode: 'copy', pattern: 'genes.tsv'
    container params.scanpy_scripts_container

    input:
        path anndata

    output:
	path 'normalised_anndata.h5ad'
	path 'matrix.mtx'
	path 'barcodes.tsv'
	path 'genes.tsv'

    script:
    """
        export PYTHONIOENCODING='utf-8'
        scanpy-normalise-data \
        --no-log-transform \
        --normalize-to '1000000.0' \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata \
        'normalised_anndata.h5ad' \
        --export-mtx ./
    """
}