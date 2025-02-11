# scxa-tertiary-workflow
Tertiary component for Single-Cell Expression Atlas workflows, focused on post-processing and advanced analyses like normalization, PCA, clustering, t-SNE, and UMAP visualizations.

## Overview

This Nextflow workflow is designed to perform downstream (tertiary) analysis on single-cell RNA sequencing (scRNA-seq) data. It takes processed data (expression matrix and metadata) as input, normalizes and scales the data, identifies variable genes, runs principal component analysis (PCA), integrates batch effects using Harmony, calculates cell neighborhoods, finds clusters, and performs visualizations like UMAP and t-SNE.

The workflow runs downstream analyses using Scanpy, leveraging the [scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts) package to run individual steps of the Scanpy workflow.

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

- Nextflow.
- SLURM.  For cluster job scheduling if using the --slurm option.
  
### Running the Workflow

The workflow can be executed for two types of scRNA-seq technologies: plate-based and droplet-based.

#### For plate-based Data:

```sh
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> [--output_path <PATH>]  [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```
#### For droplet-based data:
```sh
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> --technology droplet [--output_path <PATH>] [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```
- `--technology droplet`: Specifies that the data is droplet-based. This enables additional steps for multiplet detection (using Scrublet) and doublet removal.
- The remaining parameters are the same as for the plate-based run.

### Output

If `[--output_path <PATH>]` is not specified results will be `<EXP-ID with path>/results` dir. 
