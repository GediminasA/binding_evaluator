rule run_cadscore_eval:
    input:
        mut = work_dir + "/mutants_structure_generation/TEMPLATES/optimized_models/{pdb}={chain}={mutations}.pdb",
        wt = work_dir + "/mutants_structure_generation/TEMPLATES/optimized_models/{pdb}={chain}=nan.pdb",
        container = "containers/voronota.sif"
    output:
        work_dir + "/mutants_structure_scoring/CADscore/scores/{pdb}={chain}={mutations}.sc"
    container:
        "containers/voronota.sif"
    shell:
        "voronota-cadscore --input-target {input.wt} --input-model {input.mut} > {output}"
