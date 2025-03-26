include { SCANPY_MULTIPLET_SCRUBLET  } from '../../../modules/scanpy-scripts/scanpy_multiplet_scrublet'
include { SCANPY_PLOT_SCRUBLET } from '../../../modules/scanpy-scripts/scanpy_plot_scrublet'


process SCANPY_MULTIPLET_SCRUBLET3 {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val batch_variable

    output:
        path 'scrublet3.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
        if [ -z "$batch_variable" ]; then
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            $anndata \
            scrublet3.h5ad
        else
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            --batch-key "$batch_variable" \
            $anndata \
            scrublet3.h5ad
        fi
    """
    stub:
    """
        touch scrublet3.h5ad
    """
}

process SCANPY_MULTIPLET_SCRUBLET2 {
    container params.scanpy_scripts_container
    
    input:
        path anndata
        val batch_variable

    output:
        path 'scrublet2.h5ad'

    script:
    def args    = task.ext.args ?: ""
    """
        export PYTHONIOENCODING='utf-8'
        if [ -z "$batch_variable" ]; then
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            $anndata \
            scrublet2.h5ad
        else
            scanpy-cli multiplet scrublet \
            --input-format 'anndata' \
            --output-format 'anndata' \
            --batch-key "$batch_variable" \
            $anndata \
            scrublet3.h5ad
        fi
    """
    stub:
    """
        touch scrublet2.h5ad
    """
}


/*
 * takes collected results of different doublet finding processes and outputs majority voted results
 * requires different processes to output hd5a anndata
 */
process MAJORITY_VOTE_DOUBLET {
    container "quay.io/biocontainers/scanpy-scripts:1.1.2--pypyhdfd78af_1"
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
    def available_methods = ['scrublet1', 'scrublet2', 'scrublet'] // add new methods and tool names here
    // check if user arguments match available methods
    if (!available_methods.contains(method)) {
            error """
                Invalid doublet method: '${method}'
                Available methods: ${available_methods.join(', ')}
                """
    }
    // run all methods
    switch(method) {
        case 'scrublet1':
            return SCANPY_MULTIPLET_SCRUBLET2(adata_ch, batch_var)

        case 'scrublet2':
            return SCANPY_MULTIPLET_SCRUBLET3(adata_ch, batch_var)

        case 'scrublet':
            scrublet_result = SCANPY_MULTIPLET_SCRUBLET(adata_ch, batch_var)
            SCANPY_PLOT_SCRUBLET( scrublet_result )
            return scrublet_result 

        default:
            error "Unknown method: ${method}"
    }
}

/*
* this is heavily borrowed from nf-core/scdownstream
*/
workflow RUN_DOUBLET{

    take:
        hd_ch // channel for hd5a file
        batch_var
        methods_str
        doublet_filt
    main:
        methods = methods_str.split(',') // split methods into list
        doublet_filt_thresh = channel.value(doublet_filt)
        
        result_ch = channel.empty()
        methods.each {
            method -> result_ch = result_ch.concat(runDoubletProcess(method, hd_ch, batch_var))
        }

        // perform majority voting 
        majority_ch = MAJORITY_VOTE_DOUBLET(channel.value(params.doublet_methods), result_ch.collect(),  doublet_filt_thresh)
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
        result = majority_ch
}

workflow {
    main:
        
        hd5a = channel.fromPath(params.input_hd5a)
        RUN_DOUBLET(hd5a,
            batch_var,
            methods,
            doublet_filt_thresh)
}