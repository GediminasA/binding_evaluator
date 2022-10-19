### WILDCARDS ###
wildcard_constraints:
    forcefield = "[a-zA-Z0-9_]+"

### CONTAINER BUILDING RULES ###
rule build_openmm:
    input:
        "containers/openmm.def"
    output:
        "containers/openmm.sif"
    shell:
        "singularity build {output} {input}"

### EVALUATION RULES ###

rule evaluate_openmm:
    input:
        structure = "{directory}/{stem}_{part}_{frame}.pdb",
        container = "containers/openmm.sif"
    output:
        tsv = temp("{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_openmm_{forcefield}.tsv")
    container:
        "containers/openmm.sif"
    shell:
        """
        covid-lt/bin/pdb_openmm_minimize {input.structure} --forcefield {wildcards.forcefield}.xml --max-iterations 0 --print-forces > {output}
        """
