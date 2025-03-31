process SCANPY_FILTER_GENES {
    publishDir "${params.result_dir_path}/matrices/raw_filtered", mode: 'copy', pattern: 'matrix.mtx'
    publishDir "${params.result_dir_path}/matrices/raw_filtered", mode: 'copy', pattern: 'barcodes.tsv'
    publishDir "${params.result_dir_path}/matrices/raw_filtered", mode: 'copy', pattern: 'genes.tsv'
    container params.scanpy_scripts_container

    input:
        path anndata
        path genes

    output:
        path 'filtered_gene_anndata.h5ad'
        path 'matrix.mtx'
        path 'barcodes.tsv'
        path 'genes.tsv'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
        scanpy-filter-genes \
        --param 'g:n_cells' 3.0 1000000000.0 \
        --subset 'g:index' \
        $genes \
        --input-format 'anndata' $anndata \
        --show-obj stdout \
        --output-format anndata \
        'filtered_gene_anndata.h5ad' \
        --export-mtx ./
    """
    stub:
    """
        touch filtered_gene_anndata.h5ad
        touch matrix.mtx
        touch barcodes.tsv
        touch genes.tsv
    """
}