rule run_DSSP_part_eval:
    input:
        structure = work_dir+"/processed/{pdb}.pdb",
        container = "containers/dssp.sif"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/part/{pdb,[^=]+}={chain,[^=]+}.sc"
    container:
        "containers/dssp.sif"
    shell:
        """
        grep ^ATOM {input.structure} \
            | covid-lt-new/bin/pdb_add_header --id {wildcards.pdb} \
            | covid-lt/bin/pdb_select --chain {wildcards.chain} \
            | dssp --output-format dssp /dev/stdin > {output}
        """

rule run_DSSP_complex_eval:
    input:
        structure = work_dir+"/processed/{pdb}.pdb",
        container = "containers/dssp.sif"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{pdb,[^=]+}={chain,[^=]+}.sc"
    container:
        "containers/dssp.sif"
    shell:
        """
        grep ^ATOM {input.structure} \
            | covid-lt-new/bin/pdb_add_header --id {wildcards.pdb} \
            | dssp --output-format dssp /dev/stdin > {output}
        """

rule run_DSSP_part_sum:
    input:
        work_dir + "/mutants_structure_scoring/DSSP/scores/{type}/{pdb,[^=]+}={chain,[^=]+}.sc"
    output:
        work_dir + "/mutants_structure_scoring/DSSP/scores/{type}/{pdb,[^=]+}={chain,[^=]+}={mutations}.sum"
    shell:
        """
        join \
            <(grep -vP '\.$' {input} | awk "{{ if( \$3 == \\"{wildcards.chain}\\" ) {{ print \$2 \\" \\" substr( \$0, 36, 3 ) }} }}" | sed 's/ \+/ /g' | sort -k 1b,1) \
            <(echo {wildcards.mutations} | sed 's/+/\\n/g' | cut -c 2- | sed 's/.$//' | sort -k 1b,1) \
            | cut -d ' ' -f 2 \
            | xargs -i echo +{{}} \
            | xargs echo \
            | sed 's/^+//' \
            | bc > {output}
        """
