rule run_DSSP_part_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb,[^=]+}={chain,[^=]+}=nan.pdb",
        container = "containers/dssp.sif"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/part/{pdb,[^=]+}={chain,[^=]+}.sc"
    container:
        "containers/dssp.sif"
    shell:
        """
        covid-lt/bin/pdb_add_header --id {wildcards.pdb} {input.structure} \
            | covid-lt/bin/pdb_select --chain {wildcards.chain} \
            | dssp --output-format dssp /dev/stdin > {output}
        """

rule run_DSSP_complex_eval:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb,[^=]+}={chain,[^=]+}=nan.pdb",
        container = "containers/dssp.sif"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{pdb,[^=]+}={chain,[^=]+}.sc"
    container:
        "containers/dssp.sif"
    shell:
        """
        covid-lt/bin/pdb_add_header --id {wildcards.pdb} {input.structure} \
            | dssp --output-format dssp /dev/stdin > {output}
        """

rule run_DSSP_part_sum:
    input:
        work_dir + "/mutants_structure_scoring/DSSP/scores/part/{pdb,[^=]+}={chain,[^=]+}.sc"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/part/{pdb,[^=]+}={chain,[^=]+}={mutations}.sum"
    shell:
        """
        join -1 2 \
            <(grep -vP '\.$' {input} | awk "{{ if( \$3 == \\"{wildcards.chain}\\" ) {{ print \$0 }} }}" ) \
            <(echo {wildcards.mutations} | sed 's/+/\\n/g' | cut -c 2- | sed 's/.$//') > {output}
        """

rule run_DSSP_complex_sum:
    input:
        work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{pdb}={chain,[^=]+}.sc"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{pdb}={chain,[^=]+}={mutations}.sum"
    shell:
        """
        cp {input} {output}
        """
