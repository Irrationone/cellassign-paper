## Figure/table generation pipeline

Pipeline that generates figures/tables from the results of the analysis pipelines. 

This pipeline generates:

* Main figures
* Supplemental figures
* Supplemental tables

`config/*.yaml` contains the pipeline settings, which will be user-specific (depending on where your inputs/outputs are stored). 

This pipeline can be run with `snakemake` as:

```
snakemake -p --cores 20 --use-singularity
```

suffixed by any necessary singularity arguments (e.g. for mounting directories containing pipeline and input files, use `--singularity-args "-B local_path:container_path"`). 