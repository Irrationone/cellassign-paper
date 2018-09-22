# CellAssign simulation pipeline

Pipeline for running simulations to compare CellAssign against other methods. 

Example run:

```
toil-cwl-runner --jobStore azure:aztoil:test21 --logDebug --realTimeLogging --retryCount 0 --destBucket wasb://dest21@aztoil.blob.core.windows.net /datadrive/projects/cellassign-paper/pipelines/cellassign-sim-comparison/simulate_analysis.cwl /datadrive/projects/cellassign-paper/pipelines/cellassign-sim-comparison/config/test_azgpu.yaml
```