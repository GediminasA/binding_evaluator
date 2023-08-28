rule run_DSSP_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.pdb",
        container = "containers/dssp.sif"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/{pdb}={chain}={mutations}.sc"
    container:
        "containers/dssp.sif"
    shell:
        """
        covid-lt/bin/pdb_add_header --id {wildcards.pdb} {input.structure} \
            | dssp --output-format dssp /dev/stdin > {output}
        """
