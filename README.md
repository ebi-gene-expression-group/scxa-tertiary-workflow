# scxa-tertiary-workflow
Tertiary component for SCXA workflows

# How to run workflow for tertiary analysis 
## Prepare data
```
bash scripts/data_prep.sh <EXP-ID> [output path]
```
## Run for plate
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> [--output_path <PATH>]  [--scanpy_scripts_container <container_id>]
```
## Run for droplet
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID with path> --technology droplet [--output_path <PATH>] [--scanpy_scripts_container <container_id>]
```

If `[--output_path <PATH>]` is not specified results will be `<EXP-ID with path>/results` dir. 
