rule provean_nr:
    output:
        "data/nr/2011-08/nr.pal"
    shell:
        """
        mkdir --parents $(dirname {output})
        (
            cd $(dirname {output})
            for NUMBER in $(seq 0 5)
            do
                wget ftp://ftp.jcvi.org/data/provean/nr_Aug_2011/nr.0$NUMBER.tar.gz
                wget ftp://ftp.jcvi.org/data/provean/nr_Aug_2011/nr.0$NUMBER.tar.gz.md5
                md5sum --check nr.0$NUMBER.tar.gz.md5
                tar -xf nr.0$NUMBER.tar.gz
            done
        )
        """

# This rule is needed as containers/provean.sif does not have python3-biopython
rule run_PROVEAN_sequence:
    input:
        structure = work_dir + "/processed/{pdb}.pdb"
    output:
        work_dir + "/mutants_structure_scoring/PROVEAN/sequences/{pdb}={chain}=nan.fasta"
    shell:
        """
        covid-lt-new/bin/pdb_select --chain {wildcards.chain} {input} \
            | covid-lt-new/bin/pdb_atom2fasta > {output}
        """

rule run_PROVEAN_eval:
    input:
        sequence = work_dir + "/mutants_structure_scoring/PROVEAN/sequences/{pdb}={chain}=nan.fasta",
        container = "containers/provean.sif", # to prepare the container
        nr = "data/nr/2011-08/nr.pal" # to prepare the NR database
    output:
        work_dir + "/mutants_structure_scoring/PROVEAN/scores/{pdb}={chain}={mutations}.sc"
    log:
        work_dir + "/mutants_structure_scoring/PROVEAN/scores/{pdb}={chain}={mutations}.log"
    container:
        "containers/provean.sif"
    threads: 4
    shell:
        """
        mkdir --parents $(dirname {output})

        FASTA_FILE=$(realpath {input.sequence})
        MUT_FILE=$(mktemp)
        echo {wildcards.mutations}  | sed 's/-/del/g' | tr + '\n' > $MUT_FILE

        (cd $(dirname {input.nr}) && provean --num_threads {threads} -q $FASTA_FILE -v $MUT_FILE --psiblast psiblast --cdhit cdhit --blastdbcmd blastdbcmd -d nr) > {output} 2> {log}

        rm $MUT_FILE
        sleep 1

        """
        #echo "E144del;F145del" > $MUT_FILE
