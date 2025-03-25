process RESTORE_UNSCALED {
    container params.scanpy_scripts_container

    input:
	tuple path(anndata), path(normalise_internal_data)

    output:
	path "restore_unscaled_output_${anndata}.h5"

    script:
    """
	export PYTHONIOENCODING='utf-8'
	ln -s $anndata input.h5
	ln -s $normalise_internal_data r_source.h5
	python ${projectDir}/bin/restore_unscaled.py
	mv output.h5 'restore_unscaled_output_${anndata}.h5'
    """
}