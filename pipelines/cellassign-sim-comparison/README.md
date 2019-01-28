# Simulation analysis pipeline

Toil pipeline for running simulations to compare CellAssign against other methods. 

Example run:

```
toil-cwl-runner --jobStore {azure_job_store} --logDebug --realTimeLogging --retryCount 2 --provisioner azure --batchSystem mesos --nodeTypes Standard_A2m_v2 --maxNodes 1 --destBucket {destination_azure_bucket} simulate_analysis.cwl {config_file}
```

Docker containers containing the exact version of pipeline used for this pipeline are available at https://quay.io/repository/irrationone/toil?tab=tags. 
