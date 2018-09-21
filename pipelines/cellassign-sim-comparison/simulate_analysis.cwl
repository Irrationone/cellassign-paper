cwlVersion: v1.0
class: Workflow

requirements:
    - class: ScatterFeatureRequirement

inputs:
    script_file: File
    simulate_rmd_file: File
    cluster_rmd_file: File
    random_seeds: int[]
    num_cellss: int[]
    sim_model: string
    num_groupss: int[]
    num_batchess: int[]
    group_probs: string[]
    batch_probs: string[]
    de_faclocs: float[]
    de_facscales: float[]
    de_facnus: float[]
    de_min: float
    de_max: float
    de_probs: float[]
    down_probs: float[]
    batch_faclocs: float[]
    batch_facscales: float[]
    base_sce_params: File[]
    dimreduce_methods: string[]
    clustering_methods: string[]
    conda_env: string
    fc_percentiles: float[]
    expr_percentiles: float[]
    frac_geness: float[]
    max_geness: float[]

outputs:
    sce_cluster:
        type: File[]
        outputSource: cluster/sce_cluster

steps:
    simulate_sce:
        run: steps/create_sce.cwl
        in:
            script_file: script_file
            rmd_file: simulate_rmd_file
            random_seed: random_seeds
            num_cells: num_cellss
            sim_model: sim_model
            num_groups: num_groupss
            num_batches: num_batchess
            group_probs: group_probs
            batch_probs: batch_probs
            de_facloc: de_faclocs
            de_facscale: de_facscales
            de_facnu: de_facnus
            de_min: de_min
            de_max: de_max
            de_prob: de_probs
            down_prob: down_probs
            batch_facloc: batch_faclocs
            batch_facscale: batch_facscales
            base_sce_param: base_sce_params
        scatter: [random_seed, num_cells, num_groups, num_batches, de_facloc, de_facscale, de_facnu, de_prob, down_prob, batch_facloc, batch_facscale, base_sce_param]
        scatterMethod: flat_crossproduct

        out:
            [sce_sim]
    
    cluster:
        run: steps/cluster.cwl
        in:
            script_file: script_file
            rmd_file: cluster_rmd_file
            sce: simulate_sce/sce_sim
            dimreduce_method: dimreduce_methods
            clustering_method: clustering_methods
            conda_env: conda_env
            fc_percentile: fc_percentiles
            expr_percentile: expr_percentiles
            frac_genes: frac_geness
            max_genes: max_geness
        scatter: [sce, dimreduce_method, clustering_method, fc_percentile, expr_percentile, frac_genes, max_genes]
        scatterMethod: flat_crossproduct

        out:
            [sce_cluster]


