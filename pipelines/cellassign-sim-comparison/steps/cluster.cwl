cwlVersion: v1.0
class: CommandLineTool

hints:
    DockerRequirement:
        dockerPull: alzhang/scrna-analysis-sim-3.5

requirements:
    ResourceRequirement:
        ramMin: 8000
        tmpdirMin: 4000
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
    sce:
        type: File
        doc: SCE object to preprocess
        inputBinding:
            prefix: --sce
            position: 3
    dimreduce_method:
        type: string
        doc: Dimensionality reduction method
        inputBinding:
            prefix: --dimreduce_method
            position: 4
    clustering_method:
        type: string
        doc: Clustering method
        inputBinding:
            prefix: --clustering_method
            position: 5
    conda_env:
        type: string
        doc: Conda environment
        inputBinding:
            prefix: --conda_env
            position: 6
    fc_percentile:
        type: float
        doc: FC percentile above which to select marker genes
        inputBinding:
            prefix: --fc_percentile
            position: 7
    expr_percentile:
        type: float
        doc: Mean expr percentile above which to select marker genes
        inputBinding:
            prefix: --expr_percentile
            position: 8
    frac_genes:
        type: float
        doc: Fraction of genes above percentiles to take as marker genes
        inputBinding:
            prefix: --frac_genes
            position: 9
    max_genes:
        type: float
        doc: Maximum number of marker genes to select
        inputBinding:
            prefix: --max_genes
            position: 10
    test_proportion:
        type: float
        doc: Proportion of genes to use for testing
        inputBinding:
            prefix: --test_proportion
            position: 11


outputs:
    sce_cluster:
        type: File
        outputBinding:
            glob: "*.rds"
        secondaryFiles:
            - "^.nb.html"
            - "^_eval_measures.tsv"
            - "^_params.tsv"
            - "^_delta_compare.tsv"
    stderr_file:
        type: stderr
    stdout_file:
        type: stdout