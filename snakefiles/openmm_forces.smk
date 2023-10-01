### WILDCARDS ###
wildcard_constraints:
    ff1 = "[a-zA-Z0-9_]+",
    ff2 = "[a-zA-Z0-9_]+"

### CONTAINER BUILDING RULES ###
rule build_openmm:
    input:
        "containers/openmm.def"
    output:
        "containers/openmm.sif"
    shell:
        "apptainer build {output} {input}"

### EVALUATION RULES ###

rule evaluate_openmm_ff1:
    input:
        structure = "{directory}/{stem}_{part}_{frame}.pdb",
        container = "containers/openmm.sif"
    output:
        tsv = "{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_ff_{ff1}.tsv"
    log: 
        "{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_ff_{ff1}.log"
    container:
        "containers/openmm.sif"
    shell:
        """
        covid-lt/bin/pdb_openmm_minimize {input.structure} --forcefield {wildcards.ff1}.xml --max-iterations 0 --print-forces 2> {log}  1> {output}
        """

rule evaluate_openmm_ff2:
    input:
        structure = "{directory}/{stem}_{part}_{frame}.pdb",
        container = "containers/openmm.sif"
    output:
        tsv = "{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_ff_{ff1}-{ff2}.tsv"
    log:
        "{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_ff_{ff1}-{ff2}.log"
    container:
        "containers/openmm.sif"
    shell:
        """
        covid-lt/bin/pdb_openmm_minimize {input.structure} --forcefield {wildcards.ff1}.xml --forcefield {wildcards.ff2}.xml --max-iterations 0 --print-forces 2> {log}  1> {output}
        """

rule run_OpenMM_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/optimized_models/{pdb}={chain}={mutations}.pdb",
        groups = work_dir + "/processed_info/{pdb}_interactigGroups.tsv",
        container = "containers/openmm.sif"
    output:
        work_dir + "/mutants_structure_scoring/OpenMM/scores/{pdb}={chain}={mutations}.sc"
    container:
        "containers/openmm.sif"
    shell:
        """
        paste \
            <(sed 's/HSE/HIS/g' {input.structure} \
                | covid-lt-new/bin/pdb_openmm_minimize --forcefield charmm36.xml --forcefield implicit/gbn2.xml --print-forces --max-iterations 0 --force-unit kcal/mol --split-nonbonded-force) \
            <(sed 's/HSE/HIS/g' {input.structure} \
                | covid-lt-new/bin/pdb_select --chain $(cut -f 1 {input.groups} | sed 's/,//g') \
                | covid-lt-new/bin/pdb_openmm_minimize --forcefield charmm36.xml --forcefield implicit/gbn2.xml --print-forces --max-iterations 0 --force-unit kcal/mol --split-nonbonded-force) \
            <(sed 's/HSE/HIS/g' {input.structure} \
                | covid-lt-new/bin/pdb_select --chain $(cut -f 2 {input.groups} | sed 's/,//g') \
                | covid-lt-new/bin/pdb_openmm_minimize --forcefield charmm36.xml --forcefield implicit/gbn2.xml --print-forces --max-iterations 0 --force-unit kcal/mol --split-nonbonded-force) \
            | cut -f 1,2,4,6 > {output}
        """

rule run_OpenMM_eval_subtract:
    input:
        mut = work_dir + "/mutants_structure_scoring/OpenMM/scores/{pdb}={chain}={mutations}.sc",
        wt = work_dir + "/mutants_structure_scoring/OpenMM/scores/{pdb}={chain}=nan.sc"
    output:
        work_dir + "/mutants_structure_scoring/OpenMM/scores/{pdb}={chain}={mutations}.diff"
    shell:
        """
        paste {input.mut} {input.wt} | awk '{{ print $1 "\t" $2 - $3 - $4 - $6 + $7 + $8 }}' > {output}
        """

rule optimize_complex:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/all_models/{pdb}={chain}={mutations}.pdb",
        container = "containers/openmm.sif"
    output:
        work_dir + "/mutants_structure_generation/TEMPLATES/optimized_models/{pdb}={chain}={mutations}.pdb",
    singularity:
        "containers/openmm.sif"
    shell:
        """
        PYTHONPATH=covid-lt-new covid-lt-new/bin/pdb_renumber {input.structure} \
            | PYTHONPATH=covid-lt-new covid-lt-new/bin/pdb_resolve_alternate_locations \
            | covid-lt-new/bin/pdb_openmm_minimize --forcefield charmm36.xml --add-missing-hydrogens --constrain heavy --max-iterations 100 \
            | PYTHONPATH=covid-lt-new covid-lt-new/bin/pdb_rename_chains --source {input.structure} > {output}
        """
