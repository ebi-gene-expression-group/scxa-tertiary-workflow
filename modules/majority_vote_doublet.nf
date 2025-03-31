/*
 * takes collected results of different doublet finding processes and outputs majority voted results
 * requires different processes to output hd5a anndata
 */
process MAJORITY_VOTE_DOUBLET {
    container params.scanpy_scripts_container
    publishDir "${params.result_dir_path}/markers", mode: 'copy',  pattern: 'doublet_major*.h5ad'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    input:
        val methods 
        path h5ad
        val filter_threshold
    
    output:
        path 'doublet_major.h5ad'

    script:
    """
    majority_vote.py $methods ${h5ad.join(',')} $filter_threshold
    """
    
    stub:
    """
    touch 'doublet_major.h5ad'
    """
}