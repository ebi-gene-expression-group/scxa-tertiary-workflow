include { scanpy_plot_scrublet } from '../../../main'
include { scanpy_multiplet_scrublet } from '../../../main'
params.scanpy_scripts_container = "quay.io/biocontainers/scanpy-scripts:1.1.2--pypyhdfd78af_1"

// fake process
process doublet_method_one {
    publishDir params.result_dir_path, mode: 'copy', pattern: '.h5ad'
    input:
        path anndata
        val batch_variable

    output:
        path 'method_one.hd5a'

    script:
    """
       cp $anndata method_one.hd5a
    """
}

// fake doublet process takes h5ad output h5ad
process doublet_method_two {
    publishDir params.result_dir_path, mode: 'copy', pattern: '.h5ad'
    input:
        path anndata
        val batch_variable

    output:
        path 'method_two.hd5a'

    script:
    """
       cp $anndata method_two.hd5a
    """
}




/*
 * takes collected results of different doublet finding processes and outputs majority voted results
 * requires different processes to output hd5a anndata
 */
process MAJORITY_VOTE_DOUBLET {
    container params.scanpy_scripts_container
    publishDir params.result_dir_path, mode: 'copy', pattern: '.h5ad'
    input:
        val methods 
        path h5ad
        val filter_threshold
    
    output:
        path 'doublet_major.h5ad'

    script:
    """
    majority_vote.py $methods ${h5ad.join(',')} $filter_threshold
    echo $methods ${h5ad.join(',')} $filter_threshold >> majority_vote_params.txt
    """
}


// add new doublet finding modules in here
def runDoubletProcess(method, adata_ch, batch_var=channel.empty()) {
    def available_methods = ['one', 'two', 'scrublet'] // add new methods and tool names here
    // check if user arguments match available methods
    if (!available_methods.contains(method)) {
            error """
                Invalid doublet method: '${method}'
                Available methods: ${available_methods.join(', ')}
                """
    }
    // run all methods
    switch(method) {
        case 'one':
            return doublet_method_one(adata_ch, batch_var)

        case 'two':
            return doublet_method_two(adata_ch, batch_var)

        case 'scrublet':
            scrublet_result = scanpy_multiplet_scrublet(adata_ch, batch_var)
            scanpy_plot_scrublet( adata_ch )
            return scrublet_result 

        default:
            error "Unknown method: ${method}"
    }
}

params.batch_variable = ""
/*
* this is heavily borrowed from nf-core/scdownstream
*/
workflow RUN_DOUBLET{

    take:
        hd_ch // channel for hd5a file
        
    main:
        println params.result_dir_path
        methods = params.doublet_methods.split(',') // split methods into list
        doublet_filt_thresh = channel.value(params.doublet_filter)
        
        result_ch = channel.empty()
        methods.each {
            method -> result_ch = result_ch.concat(runDoubletProcess(method, hd_ch, params.batch_variable))
        }

        // perform majority voting 
        MAJORITY_VOTE_DOUBLET(channel.value(params.doublet_methods), result_ch.collect(),  doublet_filt_thresh)

        /* In case the pipeline is allowed to handle mutilple different datasets at once in the future

        result = methods.collectEntries {
            method ->  [(method):runDoubletProcess(methods, hd_ch, batch_val)]
        }
        result_ch = channel.empty()
        results.keySet().toSorted().each{ k -> 
          results[k].adata.collect().view()
          test = test.concat(result_ch[k].tsv)
        }
       */
        

    emit:
        result = result_ch
}

workflow {
    main:
        hd5a = channel.fromPath(params.input_hd5a)
        RUN_DOUBLET(hd5a)
}