process RUN_UMAP {
    publishDir "${params.result_dir_path}/umap", mode: 'copy', pattern: 'umap_n_neighbors_*.tsv'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    container params.scanpy_scripts_container
    
    input:
        path anndata

    output:
        path "umap_*.h5ad"
        path "umap_n_neighbors_*.tsv"

    script:
    """
	export PYTHONIOENCODING='utf-8'
	echo \$PYTHONIOENCODING
	VAR="$anndata"
	n_number="\${VAR%.h5ad}"
	echo \$n_number
	scanpy-run-umap \
            --neighbors-key "neighbors_n_\${n_number}" \
            --key-added "neighbors_n_\${n_number}" \
            --export-embedding embeddings.tsv \
            --n-components 2 \
            --min-dist 0.5 \
            --spread 1.0 \
            --alpha 1.0 \
            --gamma 1.0 \
            --negative-sample-rate 5 \
            --random-state 0 \
            --init-pos 'spectral' \
            --input-format 'anndata' \
            $anndata \
            --show-obj stdout \
            --output-format anndata \
            "umap_\${n_number}.h5ad" \
            && mv "embeddings_neighbors_n_\${n_number}.tsv" umap_n_\${n_number}.tsv

    """
}