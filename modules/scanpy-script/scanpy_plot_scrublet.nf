process SCANPY_PLOT_SCRUBLET {
    publishDir params.result_dir_path, mode: 'copy', pattern: '(scrublet.png)'
    container params.scanpy_scripts_container
    
    input:
        path anndata

    output:
        path 'scrublet.png'

    script:
    """
        export PYTHONIOENCODING='utf-8'
        scanpy-cli plot scrublet \
        --input-format "anndata" \
        --scale-hist-obs "linear" \
        --scale-hist-sim "linear" \
        $anndata \
        scrublet.png
    """
}