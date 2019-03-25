## HGSC analysis pipeline

Pipeline that runs analysis involving single cell RNA-seq data for high-grade serous ovarian cancer primary tumour samples.

This pipeline performs:

* Preprocessing and normalization
* Dimensionality reduction
* Cell cycle scoring
* CellAssign analysis
* Unsupervised clustering
* Differential expression
* Enrichment analysis

`config/*.yaml` contains the pipeline settings, which will be user-specific (depending on where your inputs/outputs are stored). 

This pipeline can be run with `snakemake` as:

```
snakemake -p --cores 20 --use-singularity
```

suffixed by any necessary singularity arguments (e.g. for mounting directories containing pipeline and input files, use `--singularity-args "-B local_path:container_path"`). 