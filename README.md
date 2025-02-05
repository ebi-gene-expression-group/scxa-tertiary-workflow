# scxa-tertiary-workflow
Tertiary component for SCXA workflows

# How to run workflow for tertiary analysis 
## Prepare data

To perform a tertiary analysis, all required datasets must be stored in a single directory, which should be specified using the --dir_path parameter.

Ensure the following files are present in the directory:

`genes_metadata.tsv` – Metadata information for genes
`genes.tsv` – List of gene identifiers
`barcodes.tsv` – Cell barcode identifiers
`matrix.mtx` – Expression matrix in Matrix Market format
`cell_metadata.tsv` – Metadata information for individual cells

## Run for plate
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> [--output_path <PATH>]  [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```
## Run for droplet
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> --technology droplet [--output_path <PATH>] [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```

If `[--output_path <PATH>]` is not specified results will be `<EXP-ID with path>/results` dir. 
