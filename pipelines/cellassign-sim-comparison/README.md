# CellAssign simulation pipeline

Pipeline for running simulations to compare CellAssign against other methods. 

Example run:

```
toil-cwl-runner --jobStore azure:aztoil:test1 --logDebug --realTimeLogging --retryCount 2 --provisioner azure --batchSystem mesos --nodeTypes Standard_A2m_v2 --maxNodes 1 --destBucket wasb://dest22@aztoil.blob.core.windows.net /root/cellassign-paper/pipelines/cellassign-sim-comparison/simulate_analysis.cwl /root/cellassign-paper/pipelines/cellassign-sim-comparison/config/test_azgpu.yaml
```