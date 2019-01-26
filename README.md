# Probabilistic cell type assignment of single-cell transcriptomic data reveals spatiotemporal microenvironment dynamics in human cancers

This repository contains the pipelines and metadata used to generate the results and figures in https://www.biorxiv.org/content/10.1101/521914v1. 


## System requirements

The pipelines under `pipelines/` can be run using [snakemake](https://snakemake.readthedocs.io/en/stable/). Most pipelines are entirely containerized, relying on public Docker images available at [my dockerhub](https://hub.docker.com/u/alzhang), and as such have no additional system requirements. These containers are quite large, as they were meant for more general purposes. 

Config files for each pipeline (`pipelines/*/config/*.yaml`) must be updated with correct data paths for pipelines to run. 

Additional dependencies:

* pymc3 3.5
* python 3.6.6
* [toil 3.18.0a1](https://github.com/Irrationone/toil)

