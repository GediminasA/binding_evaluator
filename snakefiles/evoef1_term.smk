rule run_EvoEF1_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.pdb",
        groups = work_dir + "/processed_info/{pdb}_interactigGroups.tsv",
        container = "containers/evoef1.sif"
    output:
        work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}={mutations}.sc"
    container:
        "containers/evoef1.sif"
    shell:
        "/EvoEF-master/EvoEF --command ComputeBinding --split $(cat {input.groups} | cut -f 1,2 | sed 's/,//' | sed 's/\t/,/') --pdb {input.structure} > {output}"
