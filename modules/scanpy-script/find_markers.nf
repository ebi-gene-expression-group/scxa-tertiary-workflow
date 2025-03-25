
process FIND_MARKERS {
    publishDir "${params.result_dir_path}/markers", mode: 'copy', pattern: 'markers_*.tsv'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    container params.scanpy_scripts_container

    input:
	tuple path(anndata), val(merged_group_slotname)

    output:
	path "markers_${merged_group_slotname}.h5ad"
	path "markers_*.tsv"

    script:
    """
	VAR="$merged_group_slotname"
        PREFIX="${params.slotname}_"
        echo \$VAR
        echo \$PREFIX
        if [[ "\$VAR" == *"\$PREFIX"* ]]; then
            suffix="resolution_\${VAR#\$PREFIX}"
        else
            suffix=\$VAR
        fi
        echo \$suffix

    	export PYTHONIOENCODING='utf-8'

	scanpy-find-markers \
	--save 'markers_${merged_group_slotname}.tsv' \
	--n-genes '100' \
	--groupby '${merged_group_slotname}' \
	--key-added 'markers_${merged_group_slotname}' \
	--method 'wilcoxon' \
	--use-raw  \
	--reference 'rest' \
	--filter-params 'min_in_group_fraction:0.0,max_out_group_fraction:1.0,min_fold_change:1.0'  \
	--input-format 'anndata' \
	$anndata  \
	--show-obj stdout \
	--output-format anndata \
	"markers_${merged_group_slotname}.h5ad" 

	command_exitcode=\$?
	echo "Command exit code: \$command_exitcode"
	
	if [ "\$command_exitcode" -eq 0 ]; then
	    if [ "${merged_group_slotname}" != "\${suffix}" ]; then
	        mv "markers_${merged_group_slotname}.tsv" "markers_\${suffix}.tsv"
	        echo "Renamed markers file to markers_\${suffix}.tsv"
	    else
	        echo "${merged_group_slotname} and \${suffix} are the same, renaming file not required."
	    fi
	else
	    echo "scanpy-find-markers failed with exit code \$command_exitcode" >&2
	    exit \$command_exitcode
	fi
    """
}