## Test datasets for Single-Cell Expression Atlas workflow

The workflow can be executed for two types of scRNA-seq technologies: plate-based and droplet-based.

For plate-based data (test-plate) (E-MTAB-9801):

```
nextflow run main.nf --slurm -resume --dir_path test-plate test-plate-out --celltype_field authors_cell_type_-_ontology_labels
```

For droplet-based data (test-droplet) (E-GEOD-130148):

```
nextflow run main.nf --slurm -resume --dir_path test-droplet test-droplet-out --technology droplet --celltype_field authors_cell_type_-_ontology_labels
```
here, `--celltype_field` is optional 
