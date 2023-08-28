rule run_EvoEF1_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.pdb",
        groups = work_dir + "/processed_info/{pdb}_interactigGroups.tsv"
    output:
        work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}={mutations}.sc"
    shell:
        "EvoEF --command ComputeBinding --split $(cat {input.groups} | cut -f 1,2 | sed 's/,//' | sed 's/\+/,/') --pdb {input.structure} > {output}"
