## Test dataset 

The workflow can be executed for two types of scRNA-seq technologies: plate-based and droplet-based.

For plate-based data (test-1):

```
nextflow run main.nf --slurm -resume --dir_path test-1 test-1-out --celltype_field authors_cell_type_-_ontology_labels
```

For droplet-based data (test-2):

```
nextflow run main.nf --slurm -resume --dir_path test-1 test-1-out --technology droplet --celltype_field authors_cell_type_-_ontology_labels
```

