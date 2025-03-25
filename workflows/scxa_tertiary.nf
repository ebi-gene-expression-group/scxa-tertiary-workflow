/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { COLUMN_REARRANGE_1            } from "${projectDir}/modules/scanpy-script/column_rearrange_1.nf"
include { COLUMN_REARRANGE_2            } from "${projectDir}/modules/scanpy-script/column_rearrange_2.nf"
include { MERGEGENEFILES                } from "${projectDir}/modules/scanpy-script/merge_gene_files.nf"
include { SCANPY_READ_10X               } from "${projectDir}/modules/scanpy-script/scanpy_read_10x.nf"
include { SCANPY_MULTIPLET_SCRUBLET     } from "${projectDir}/modules/scanpy-script/scanpy_multiplet_scrublet.nf"
include { SCANPY_PLOT_SCRUBLET          } from "${projectDir}/modules/scanpy-script/scanpy_plot_scrublet.nf"
include { SCANPY_FILTER_CELLS           } from "${projectDir}/modules/scanpy-script/scanpy_filter_cells.nf"
include { SCANPY_FILTER_GENES           } from "${projectDir}/modules/scanpy-script/scanpy_filter_genes.nf"
include { NORMALISE_DATA                } from "${projectDir}/modules/scanpy-script/normalise_data.nf"
include { NORMALISE_INTERNAL_DATA       } from "${projectDir}/modules/scanpy-script/normalise_internal_data.nf"
include { FIND_VARIABLE_GENES           } from "${projectDir}/modules/scanpy-script/find_variable_genes.nf"
include { SCALE_DATA                    } from "${projectDir}/modules/scanpy-script/scale_data.nf"
include { RUN_PCA                       } from "${projectDir}/modules/scanpy-script/run_pca.nf"
include { HARMONY_BATCH                 } from "${projectDir}/modules/scanpy-script/harmony_batch.nf"
include { NEIGHBORS                     } from "${projectDir}/modules/scanpy-script/neighbors.nf"
include { NEIGHBORS_FOR_UMAP            } from "${projectDir}/modules/scanpy-script/neighbors_for_umap.nf"
include { RUN_TSNE                      } from "${projectDir}/modules/scanpy-script/run_tsne.nf"
include { RUN_UMAP                      } from "${projectDir}/modules/scanpy-script/run_umap.nf"
include { FIND_CLUSTERS                 } from "${projectDir}/modules/scanpy-script/find_clusters.nf"
include { RESTORE_UNSCALED              } from "${projectDir}/modules/scanpy-script/restore_unscaled.nf"
include { FIND_MARKERS                  } from "${projectDir}/modules/scanpy-script/find_markers.nf"
include { MAKE_PROJECT_FILE             } from "${projectDir}/modules/scanpy-script/make_project_file.nf"

workflow SCXA_TERTIARY {

    // Create input channel (single file via CLI parameter)
    genemeta                    = Channel.fromPath("${params.dir_path}/genes_metadata.tsv")
    genes                       = Channel.fromPath("${params.dir_path}/genes.tsv")
    barcodes                    = Channel.fromPath("${params.dir_path}/barcodes.tsv")
    matrix                      = Channel.fromPath("${params.dir_path}/matrix.mtx")
    cellmeta                    = Channel.fromPath("${params.dir_path}/cell_metadata.tsv")
    neighbors_ch                = channel.fromList(params.neighbor_values)
    perplexity_ch               = channel.fromList(params.perplexity_values)
    resolution_ch               = channel.fromList(params.resolution_values)
    merged_group_slotname_ch    = Channel.fromList(params.merged_group_slotname)

    COLUMN_REARRANGE_1(
        genemeta, 
        "gene_id"
    )
    COLUMN_REARRANGE_2(
        genemeta, 
        "gene_id", 
        "gene_name"
    )
    MERGEGENEFILES(
        genes,
        COLUMN_REARRANGE_2.out
    )
    SCANPY_READ_10X(
        matrix,
        MERGEGENEFILES.out,
        barcodes,
        cellmeta,
        genemeta
    )

    if ( params.technology == "droplet" ) {
        SCRUBLET_ch = SCANPY_MULTIPLET_SCRUBLET(
            SCANPY_READ_10X.out,
            params.batch_variable
        )
        SCANPY_PLOT_SCRUBLET(
            SCRUBLET_ch
        )
        SCANPY_FILTER_CELLS(
            SCRUBLET_ch,
            "--category predicted_doublet False"
        )
    }
    else {
        SCANPY_FILTER_CELLS(
            SCANPY_READ_10X.out,
            ""
        )
    }

    SCANPY_FILTER_GENES(
        SCANPY_FILTER_CELLS.out,
        COLUMN_REARRANGE_1.out[0]
    )
    NORMALISE_DATA(
        SCANPY_FILTER_GENES.out[0]
    )
    NORMALISE_INTERNAL_DATA(
        SCANPY_FILTER_GENES.out[0]
    )
    FIND_VARIABLE_GENES(
        NORMALISE_INTERNAL_DATA.out,
        params.batch_variable
    )

    if ( params.technology == "droplet" ) {
        SCALE_DATA(
            FIND_VARIABLE_GENES.out
        )
        RUN_PCA(
            SCALE_DATA.out
        )
    }
    else {
        RUN_PCA(
            FIND_VARIABLE_GENES.out
        )
    }

    HARMONY_BATCH(
        RUN_PCA.out,
        params.batch_variable
    )
    NEIGHBORS(
        HARMONY_BATCH.out,
        params.representation
    )
    NEIGHBORS_FOR_UMAP(
        HARMONY_BATCH.out.combine(neighbors_ch),
        params.representation
    )
    TNSEs_ch = RUN_TSNE(
        HARMONY_BATCH.out.combine(perplexity_ch),
        params.representation
    )[0]

    UMAPs_ch = RUN_UMAP(
        NEIGHBORS_FOR_UMAP.out.flatten()
    )[0]

    FIND_CLUSTERS(
        NEIGHBORS.out.combine(resolution_ch)
    )

    // Combine the outputs of FIND_CLUSTERS and NEIGHBORS processes
    combined_outputs = FIND_CLUSTERS.out[0].mix(NEIGHBORS.out)

    if ( params.technology == "droplet" ) {
        RESTORE_UNSCALED(
            combined_outputs.combine(NORMALISE_INTERNAL_DATA.out)
        )
        restore_unscaled_files = RESTORE_UNSCALED.out.map { file ->
            // Extract the sample number from the file name
            def sampleNumber = file.baseName.replaceFirst('restore_unscaled_output_', '').replaceFirst('clusters', params.slotname).replaceFirst('neighbors',params.celltype_field).replaceFirst('.h5ad','')
            [file, sampleNumber] // Create a tuple with sample number and file
        }
        FIND_MARKERS(
            restore_unscaled_files
        )
    }
    else {
        processed_files = combined_outputs.map { file ->
            // Extract the sample number from the file name
            def sampleNumber = file.baseName.replaceFirst('clusters', params.slotname).replaceFirst('neighbors',params.celltype_field)
            [file, sampleNumber] // Create a tuple with sample number and file
        }
        FIND_MARKERS(
            processed_files
        )
    }
    MAKE_PROJECT_FILE(
        NEIGHBORS.out,
        SCANPY_READ_10X.out,
        SCANPY_FILTER_GENES.out[0],
        NORMALISE_DATA.out[0],
        FIND_MARKERS.out[0].collect(),
        TNSEs_ch.mix(UMAPs_ch).collect()
    )
}