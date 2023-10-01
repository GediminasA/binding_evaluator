rule run_EvoEF1_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/optimized_models/{pdb}={chain}={mutations}.pdb",
        groups = work_dir + "/processed_info/{pdb}_interactigGroups.tsv",
        container = "containers/evoef1.sif"
    output:
        work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}={mutations}.sc"
    container:
        "containers/evoef1.sif"
    shell:
        "/EvoEF-master/EvoEF --command ComputeBinding --split $(cat {input.groups} | cut -f 1,2 | sed 's/,//' | sed 's/\t/,/') --pdb {input.structure} > {output}"

rule run_EvoEF1_eval_subtract:
    input:
        mut = work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}={mutations}.sc",
        wt = work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}=nan.sc"
    output:
        work_dir + "/mutants_structure_scoring/EvoEF1/scores/{pdb}={chain}={mutations}.diff"
    shell:
        """
        grep -P '^Total\s+=' --no-filename {input.mut} {input.wt} \
            | xargs echo \
            | awk '{{ print $3 - $6 }}' > {output}
        """
