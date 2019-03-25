## Model fitting analysis pipeline

Pipeline that fits the modified splatter model to logFC values from comparing different cell types. 

This pipeline:

* Fits models with `pymc3`
* Diagnoses model fits
* Generates plots

Unlike the other pipelines, this pipeline is not Dockerized and requires installation of `pymc3` (tested on v3.5). 

In situations where disk speed appears to limit runtime, switching the `theano` compliation directory to one on a RAM drive may help (e.g. `/dev/shm`). 

`config/*.yaml` contains the pipeline settings, which will be user-specific (depending on where your inputs/outputs are stored). 

This pipeline can be run with `snakemake` as:

```
snakemake -p --cores 20 --use-singularity
```

suffixed by any necessary singularity arguments (e.g. for mounting directories containing pipeline and input files, use `--singularity-args "-B local_path:container_path"`). 