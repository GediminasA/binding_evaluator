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

rule run_PROVEAN_eval:
    input:
        sequence = work_dir + "/processed_info/{stem,[^_]+}.fasta",
        container = "containers/provean.sif", # to prepare the container
        nr = "data/nr/2011-08/nr.pal" # to prepare the NR database
    output:
        work_dir + "/mutants_structure_scoring/PROVEAN/scores/{pdb}={chain}={mutations}.sc"
    container:
        "containers/provean.sif"
    shell:
        """
        mkdir --parents $(dirname {output})

        FASTA_FILE=$(realpath {input.sequence})
        MUT_FILE=$(mktemp)

        echo {wildcards.mutations} > $MUT_FILE

        (cd $(dirname {input.nr}) && provean -q $FASTA_FILE -v $MUT_FILE --psiblast psiblast --cdhit cdhit --blastdbcmd blastdbcmd -d nr) > {output}

        rm $MUT_FILE
        """
