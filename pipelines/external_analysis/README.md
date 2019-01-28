## Koh et al. analysis pipeline

Pipeline that runs analysis involving Koh et al. (2016) data. 

This pipeline performs:

* Dimensionality reduction (using preprocessed data from Duo et al. 2018)
* CellAssign analysis
* Unsupervised clustering

`config/*.yaml` contains the pipeline settings, which will be user-specific (depending on where your inputs/outputs are stored). 

This pipeline can be run with `snakemake` as:

```
snakemake -p --cores 20 --use-singularity
```

suffixed by any necessary singularity arguments (e.g. for mounting directories containing pipeline and input files, use `--singularity-args "-B local_path:container_path"`). 
