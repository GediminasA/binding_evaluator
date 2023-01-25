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
    singularity:
        "containers/openmm.sif"
    shell:
        """
        covid-lt/bin/pdb_openmm_minimize {input.structure} --forcefield {wildcards.ff1}.xml --max-iterations 0 --print-forces > {output}
        """

rule evaluate_openmm_ff2:
    input:
        structure = "{directory}/{stem}_{part}_{frame}.pdb",
        container = "containers/openmm.sif"
    output:
        tsv = "{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_ff_{ff1}-{ff2}.tsv"
    singularity:
        "containers/openmm.sif"
    shell:
        """
        covid-lt/bin/pdb_openmm_minimize {input.structure} --forcefield {wildcards.ff1}.xml --forcefield {wildcards.ff2}.xml --max-iterations 0 --print-forces #> {output}
        """
