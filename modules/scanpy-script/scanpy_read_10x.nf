process SCANPY_READ_10X {
    container params.scanpy_scripts_container
    
    input:
        path matrix
        path genes
        path barcodes
        path cellmeta
        path genemeta

    output:
        path 'anndata.h5ad'

    script:
    """
        #ln -s $matrix matrix.mtx
        ln -s $genes genes.tsv
        #ln -s $barcodes barcodes.tsv

        export PYTHONIOENCODING='utf-8'
        
        scanpy-read-10x --input-10x-mtx ./ \
        --var-names 'gene_ids' \
        --extra-obs $cellmeta \
        --extra-var $genemeta \
        --show-obj stdout \
        --output-format anndata \
        'anndata.h5ad'
    """
}