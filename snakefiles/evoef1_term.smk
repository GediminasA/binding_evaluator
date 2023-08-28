rule run_EvoEF1_eval:
    input:
        structure = ""
    output:
        work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}={mutations}.sc"
    shell:
        "EvoEF --command ComputeBinding --split $(echo $BASE | cut -d _ -f 3-4 | sed 's/_/,/g') --pdb {input.structure} > {output}"
