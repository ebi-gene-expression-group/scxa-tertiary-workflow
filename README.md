# scxa-tertiary-workflow
Tertiary component for SCXA workflows

# How to run workflow for tertiary analysis 
## Prepare data
```
bash scripts/data_prep.sh <EXP-ID>
```
## Run for plate
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID>
```
## Run for droplet
```
nextflow run main.nf --slurm -resume --dir_path <EXP-ID> --technology droplet
```
