## Test datasets for Single-Cell Expression Atlas workflow

The workflow can be executed for two types of scRNA-seq technologies: plate-based and droplet-based.

1. For plate-based data (test-plate). More details on SC Expression Atlas [E-GEOD-9801](https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-9801):

```
nextflow run main.nf \
--slurm -resume \
--dir_path test-plate \
--output_path test-plate-out \
--celltype_field authors_cell_type_-_ontology_labels
```
here, `--technology` is not specified because it set to plate by default.

2. For droplet-based data (test-droplet). More details on SC Expression Atlas [E-GEOD-130148](https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-130148):

```
nextflow run main.nf \
--slurm -resume \
--dir_path test-droplet \
--output_path test-droplet-out \
--technology droplet \
--celltype_field authors_cell_type_-_ontology_labels
```
here, `--celltype_field` is optional 
