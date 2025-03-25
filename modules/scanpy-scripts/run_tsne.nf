process RUN_TSNE {
    publishDir "${params.result_dir_path}/tsne", mode: 'copy', pattern: 'tsne_perplexity_*\\.tsv'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    
    container params.scanpy_scripts_container
    
    input:
        tuple path(anndata), val(perplexity_values)
        val representation

    output:
        path "tsne_${perplexity_values}.h5ad"
	    path "tsne_perplexity_${perplexity_values}.tsv"

    script:
    def args    = task.ext.args ?: ""
    """
            export PYTHONIOENCODING='utf-8'
            scanpy-run-tsne \
            --use-rep $representation \
            --export-embedding embeddings.tsv \
            --perplexity $perplexity_values \
            --key-added 'perplexity_$perplexity_values' \
            --early-exaggeration '12.0' \
            --learning-rate '400.0' \
            --no-fast-tsne \
            --random-state 1234  \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            'tsne_${perplexity_values}.h5ad' \
            && mv 'embeddings_perplexity_${perplexity_values}.tsv' 'tsne_perplexity_${perplexity_values}.tsv'
    """
    stub:
    """
        touch "tsne_${perplexity_values}.h5ad"
        touch "tsne_perplexity_${perplexity_values}.tsv"
    """
}