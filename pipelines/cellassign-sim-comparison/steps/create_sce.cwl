cwlVersion: v1.0
class: CommandLineTool

hints:
    DockerRequirement:
        dockerPull: alzhang/scrna-analysis-sim-3.5

requirements:
    ResourceRequirement:
        ramMin: 16000
        tmpdirMin: 14000
        outdirMin: 4000
        coresMin: 1

baseCommand: Rscript

inputs:
    script_file:
        type: File
        doc: Wrapper script file to run Rmarkdown file from
        inputBinding:
            position: 1
    rmd_file:
        type: File
        doc: Rmarkdown file to generate report from
        inputBinding:
            prefix: --input_file
            position: 2
    random_seed:
        type: int
        doc: Random seed for simulation
        inputBinding:
            prefix: --random_seed
            position: 3
    num_cells:
        type: int
        doc: Number of cells to simulate
        inputBinding:
            prefix: --num_cells
            position: 4
    sim_model:
        type: string
        doc: Model to simulate under. Use 'splat'. Compatibility uncertain otherwise. 
        inputBinding:
            prefix: --sim_model
            position: 5
    num_groups:
        type: int
        doc: Number of groups to simulate
        inputBinding:
            prefix: --num_groups
            position: 6
    num_batches:
        type: int
        doc: Number of batches to simulate
        inputBinding:
            prefix: --num_batches
            position: 7
    group_probs: 
        type: string
        doc: Group probabilities
        inputBinding:
            prefix: --group_probs
            position: 8
    batch_probs:
        type: string
        doc: Batch probabilities
        inputBinding:
            prefix: --batch_probs
            position: 9
    de_facloc:
        type: float
        doc: Location parameter for DE
        inputBinding:
            prefix: --de_facloc
            position: 10
    de_facscale:
        type: float
        doc: Scale parameter for DE
        inputBinding:
            prefix: --de_facscale
            position: 11
    de_nu:
        type: float
        doc: Dispersion parameter for DE
        inputBinding:
            prefix: --de_nu
            position: 12
    de_prob:
        type: float
        doc: Probability of a gene being differentially expressed
        inputBinding:
            prefix: --de_prob
            position: 13
    down_prob:
        type: float
        doc: Probability of a gene being downregulated
        inputBinding:
            prefix: --down_prob
            position: 14
    de_min:
        type: float
        doc: Min FC of a gene (prior to reciprocal trans. for down_prob)
        inputBinding:
            prefix: --de_min
            position: 15
    de_max:
        type: float
        doc: Max FC of a gene (prior to reciprocal trans. for down_prob)
        inputBinding:
            prefix: --de_max
            position: 16
    batch_facloc:
        type: float
        doc: Location parameter for batch
        inputBinding:
            prefix: --batch_facloc
            position: 17
    batch_facscale:
        type: float
        doc: Scale parameter for batch
        inputBinding:
            prefix: --batch_facscale
            position: 18
    base_sce_param:
        type: File
        doc: Base SCE params object
        inputBinding:
            prefix: --base_sce_param 
            position: 19

outputs:
    sce_sim:
        type: File
        outputBinding:
            glob: "*.rds"
        secondaryFiles:
            - "^.html"
    
    stderr_file:
        type: stderr
    stdout_file:
        type: stdout