# scxa-tertiary-workflow
Tertiary component for SCXA workflows

# How to run workflow for tertiary analysis 
## Prepare data

Datasets for the tertiary analysis must be in a single directory to be used in `--dir_path`

Following files are required:
1. genes_metadata.tsv
2. genes.tsv
3. barcodes.tsv
4. matrix.mtx
5. cell_metadata.tsv

## Run for plate
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> [--output_path <PATH>]  [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```
## Run for droplet
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> --technology droplet [--output_path <PATH>] [--scanpy_scripts_container <container_id>] [--celltype_field <celltype_field>]
```

If `[--output_path <PATH>]` is not specified results will be `<EXP-ID with path>/results` dir. 
