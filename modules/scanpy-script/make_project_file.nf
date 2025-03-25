process MAKE_PROJECT_FILE {
    publishDir params.result_dir_path, mode: 'copy'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    container params.scanpy_scripts_container

    input:
        path neighbors
        path scanpy_read_10x
        path filter_genes
        path normalise_data
        path find_markers
        path TNSEs_mix_UMAPs
    output:
        path "project.h5ad"
    script:
    """
        export PYTHONIOENCODING='utf-8'
        ln -s $neighbors input.h5
        ln -s $scanpy_read_10x r_source.h5
        ln -s '$filter_genes' x_source_0.h5
        ln -s '$normalise_data' x_source_1.h5
        count=0
        for i in $find_markers
        do
                ln -sf "\${i}" obs_source_\${count}.h5
                ln -sf "\${i}" uns_source_\${count}.h5
                count=\$((count + 1))
                echo "\${count}"
        done
        count=0
        for i in $TNSEs_mix_UMAPs
        do
                ln -sf "\${i}" embedding_source_\${count}.h5
                count=\$((count + 1))
                echo "\${count}"
        done
        python "${projectDir}/bin/final_project.py"
        mv output.h5 project.h5ad
    """
}