## Runtime benchmarking pipeline

Pipeline that benchmarks runtime of cellassign on multi-group data simulated with `splatter`. 

`config/*.yaml` contains the pipeline settings, which will be user-specific (depending on where your inputs/outputs are stored). 

This pipeline can be run with `snakemake` as:

```
snakemake -p --cores 20 --use-singularity
```

suffixed by any necessary singularity arguments (e.g. for mounting directories containing pipeline and input files, use `--singularity-args "-B local_path:container_path"`). 
