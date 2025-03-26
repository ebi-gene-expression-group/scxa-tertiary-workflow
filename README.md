# scxa-tertiary-workflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)


## Introduction
Tertiary component for [Single-Cell Expression Atlas](http://www.ebi.ac.uk/gxa/sc/) workflows, focused on post-processing and advanced analyses like normalization, PCA, clustering, t-SNE, and UMAP visualizations.

## Overview

This Nextflow workflow is designed to perform analysis downstream of the quantification of expression counts from single-cell RNA sequencing (scRNA-seq) raw data. This tertiary analysis takes processed data (expression matrix and metadata) as input, normalizes and scales the data, identifies variable genes, runs principal component analysis (PCA), integrates batch effects using Harmony, calculates cell neighborhoods, finds clusters, and performs visualizations like UMAP and t-SNE.  It also finds markers for cell groupings.

The workflow runs these analyses using Scanpy, leveraging the [scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts) package to run individual steps of the Scanpy workflow.

## How to run the workflow

### Prepare the data

To perform a tertiary analysis, all required datasets must be stored in a single directory, which should be specified using the `--dir_path` parameter.

Ensure the following files are present in the directory:

- `genes_metadata.tsv` – Metadata information for genes
- `genes.tsv` – List of gene identifiers
- `barcodes.tsv` – Cell barcode identifiers
- `matrix.mtx` – Expression matrix in Matrix Market format
- `cell_metadata.tsv` – Metadata information for individual cells

### Requirements

- Nextflow
- Singularity or Docker
  
### High-performance computing

This workflow can be run on High-performance computing.

- SLURM.  For SLURM job scheduling - use the `--slurm` option.
  
### Running the workflow

The workflow can be executed for two types of scRNA-seq technologies: plate-based and droplet-based.

#### For plate-based data:

```sh
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> [--output_path <PATH>]  [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```
#### For droplet-based data:
```sh
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> --technology droplet [--output_path <PATH>] [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```
- `--technology droplet`: Specifies that the data is droplet-based. This enables additional steps for multiplet detection (using Scrublet) and doublet removal.
- The remaining parameters are the same as for the plate-based run.

#### Running the workflow for SCEA

If running for Single-cell Expression Atlas, specify the Atlas-specific profile by adding `-profile atlas` to the Nextflow command.


### Output

If `[--output_path <PATH>]` is not specified results will be `<EXP-ID with path>/results` dir. 


## Credits

ebi-gene-expression-group/scxa-tertiary-workflow was originally written by [Anil Thanki](https://github.com/anilthanki), [Iris Yu](https://github.com/irisdianauy) and [Pedro Madrigal](https://github.com/pmb59), based on a previous Galaxy workflow developed by Pablo Moreno and Jonathan Manning.

We thank the following people and teams for their extensive assistance in the development of this pipeline:

- [Khoi Bui Dinh](https://github.com/imkhoibui)
- [Jose Marugan](https://github.com/jcmarca)
- [Arsenios Chatzigeorgiou](https://github.com/arschat)

## Citations

For now, if you use the workflow for your analysis please cite it using the following doi: [10.1093/nar/gkad1021](https://doi.org/10.1093/nar/gkad1021)

## Acknowledging nf-core

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:



