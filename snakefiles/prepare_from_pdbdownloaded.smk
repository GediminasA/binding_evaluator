import os
### CONTAINER BUILDING RULES ###

rule build_promod:
    input:
        "containers/promod.def"
    output:
        "containers/promod.sif"
    shell:
        "apptainer build {output} {input}"




### EVALUATION RULES ###






# returns chains remaining  without a main chain
# if the expected main chain doesn't exits then the largest chain is taken as the main chain

pdbproc_dir = work_dir + "/pdb_proc"

rule download_pdb:
    output:
        pdbproc_dir+"/pristine/{stem}.pdb"
        #config["structures_folder"]+"/{stem}.pdb"
    log:
        pdbproc_dir+"/pristine/{stem}.log"
    params:
        stem = "{stem}"
    params:
        localdb = config["structures_folder"]
    run:
        logf = open(f"{log}","w")
        local = config["structures_folder"]+"/"+wildcards.stem+".pdb"
        if os.path.exists(local):
            logf.write(f"copying file {local} to {output}\n")
            shell("cp {local} {output}")
        else:
            logf.write(f"downloading {params.stem} to {output}\n")
            shell("pdb_fetch {params.stem} > {output}")
        logf.close()

rule download_swiss_prot_db:
    params:
        link = "https://ftp.ncbi.nlm.nih.gov/blast/db/v5/swissprot.tar.gz",
        diro = "data/blast/swiss_prot",
        archive = "swiss.tar.gz",
    output:
        swissdb = "data/blast/swiss_prot/swissprot.pin"
    shell:
        """
            mkdir -p {params.diro}
            wget {params.link} -O {params.diro}/{params.archive}
            cd {params.diro}
            tar xvzf {params.archive}
        """

rule extract_seqs_initial:
    input:
        pdbproc_dir+"/pristine/{stem}.pdb"
    output:
        fasta = work_dir+"/initial_cleanup/{stem,[^_]+}.fasta",
        pdb = work_dir+"/initial_cleanup/{stem,[^_]+}_woHetinSeqres.pdb"
    log:
         work_dir+"/initial_cleanup/{stem,[^_]+}_woHetinSeqres.log"
    params:
        path = "covid-lt/"
    # notebook:
    #     "notebooks/extract_pdbsequence.py.ipynb"
    script: 
        "notebooks/extract_pdbsequence.py.py"

rule extract_seqs_pristine:
    input:
        pdbproc_dir+"/pristine/{stem}.pdb"
    output:
        pdbproc_dir+"/pristine/{stem,[^_]+}_info.txt"
    shell: 
      """
        pdb_wc  {input} > {output}
      """
#to fix seqres entries with Z letters and alike
rule search_against_swiss_prot:
    input:
        db = "data/blast/swiss_prot/swissprot.pin",
        match = work_dir+"/initial_cleanup/{stem}.fasta"
    output:
        work_dir+"/initial_cleanup/{stem,[^_]+}_swissmatch.txt"
    params:
        db = "data/blast/swiss_prot/swissprot"
    threads: 2
    shell:
        """
            blastp -query {input.match} -db {params.db} -num_threads  {threads}\
            -out {output} -outfmt "6 evalue pident qstart qend qseqid sseqid qseq sseq"

        """

rule fix_ambiquous_positions:
    input:
        swiss = work_dir+"/initial_cleanup/{stem}_swissmatch.txt",
        pdb = work_dir+"/initial_cleanup/{stem,[^_]+}_woHetinSeqres.pdb"
    output:
        work_dir+"/initial_cleanup/{stem,[^_]+}_noambi.pdb"
    log:
        work_dir+"/initial_cleanup/{stem,[^_]+}_noambi.log"
    params:
        path = "covid-lt/"
    # notebook:
    #    "notebooks/fix_ambi.py.ipynb"
    script:
       "notebooks/fix_ambi.py.py"

rule t4t:
    input:
        expand(
            #pdbproc_dir+"/pristine/{stem}_info.txt",
            #work_dir+"/processed/{stem}.pdb",
            work_dir+"/static/splits/{stem}_0_part1.pdb",
            #work_dir+"/processed_info/{stem}_interactigGroups.tsv",
            #work_dir+"/initial_cleanup/{stem}_noambi.pdb",
            #work_dir+"/initial_cleanup/{stem}_swissmatch.txt",
            stem=pdb_stems)
            #stem=["4CPA"])




rule renumber_against_seqres:
    input:
        work_dir+"/initial_cleanup/{stem}_noambi.pdb",
        #pdbproc_dir  + "/pristine/{stem}.pdb"
    output:
        pdbproc_dir + "/process/{stem,[^_]+}_seqresMatched.pdb"
    log:
        pdbproc_dir + "/process/{stem}_seqresMatched.log"
    shell:
        """
        PYTHONPATH=covid-lt covid-lt/bin/pdb_align {input} 1> {output} 2> {log}
        """ 


rule fix_with_promod:
    input:
        "{stem}.pdb",
        "containers/promod.sif"
    output:
        "{stem}_promod.pdb"
    log:
        "{stem}_promod.log"
    container: "containers/promod.sif"
    threads: 4
    shell:
        """
        export OPENMM_CPU_THREADS={threads}
        PYTHONPATH=covid-lt covid-lt/bin/promod-fix-pdb --trim {input[0]} 1> {output} 2> {log}
        """ 


rule pdb2pqr:
    input:
        "{stem}.pdb"
    output:
        pqr = "{stem}_pdb2pqr.pqr",
        pdb = "{stem}_pdb2pqr.pdb"
    params:
        pdb = "{stem}_pdb2pqrWOseqres.pdb",
    log:
        "{stem}_pdb2pqr.log",
    shell:
        """
        grep -e "^SEQRES" {input} > {output.pdb}
        pdb2pqr30 --drop-water --include-header  {input} {output.pqr}  --pdb-output {params.pdb} &> {log}
        cat {params.pdb} >> {output.pdb}
        rm {params.pdb}
        """

rule removeHOH:
    input:
        "{stem}.pdb"
    output:
        pdb = "{stem}_woHOH.pdb",
    shell:
        """
            pdb_delresname -HOH {input} > {output}
        """

rule removeHET:
    input:
        "{stem}.pdb"
    output:
        pdb = "{stem}_woHET.pdb",
    shell:
        """
            pdb_delhetatm {input} > {output}
        """

rule removeH:
    input:
        "{stem}.pdb"
    output:
        pdb = "{stem}_woH.pdb",
    shell:
        """
            pdb_delelem -H {input} > {output}
        """
 
rule pdb_fix:
    input:
        "{stem}.pdb"
    output:
        pdb = "{stem}_pdbfix.pdb"
    threads: 24
    params:
        pdb = "{stem}_pdbfix_tmp.pdb"
    shell:
        """
        grep -e "^SEQRES" {input} > {output.pdb}
        export OPENMM_CPU_THREADS={threads}
        pdbfixer {input} --replace-nonstandard --add-atoms=heavy --add-residues --output={params.pdb}
        cat {params.pdb} >> {output.pdb}
        rm {params.pdb}
        """

rule copy_alreadyprepared4fix:
        input:
            config["preprocessed_structures"]+"/{stem}.pdb"
        output:
            work_dir+"/process_provided/{stem,[^_]+}.pdb"
        shell:
            "cp {input} {output}"

rule copy_alreadyprepared:
        input:
            work_dir+"/process_provided/{stem}.pdb"
            #work_dir+"/process_provided/{stem}_pdb2pqr.pdb"
        output:
            work_dir+"/processed/{stem,[^_]+}.pdb"
        shell:
            "cp {input} {output}"

rule copy_prepared:
        input:
            pdbproc_dir + "/process/{stem}_seqresMatched_woHOH_pdbfix_pdb2pqr_promod.pdb"
            #pdbproc_dir + "/process/{stem}_seqresMatched_woHOH_pdbfix_promod.pdb"
        output:
            work_dir+"/processed/{stem,[^_]+}.pdb"
        shell:
            "cp {input} {output}"

ruleorder:  copy_alreadyprepared > copy_prepared

#extract fasta files and run hhblitz to identify antibodies

rule extract_seqs:
    input:
        work_dir+"/processed/{stem}_woHET.pdb"
    output:
        work_dir+"/processed_info/{stem,[^_]+}.fasta"
    shell:
        "pdb_tofasta -multi {input} > {output} "

rule extract_seqs_by_chain:
    input:
        structure = pdbproc_dir + "/pristine/{stem,[^_]+}.pdb",
        template = config["mutants_templates"],
        container = "containers/muscle.sif"
    output:
        work_dir+"/processed_info/{stem}_chain_{chain}.fasta"
    container:
        "containers/muscle.sif"
    shell:
        """
        (
            covid-lt-new/bin/pdb_add_header --id {wildcards.stem} {input.structure} \
                | covid-lt-new/bin/pdb_atom2fasta \
                | covid-lt-new/bin/fasta_select --id {wildcards.stem}:{wildcards.chain}
            echo '>template'
            tail -n +2 {input.template}
        ) \
            | muscle \
            | covid-lt-new/bin/fasta_select --id {wildcards.stem}:{wildcards.chain} > {output}
        """

rule search4antibodies:
    input:
        work_dir+"/processed_info/{stem}.fasta"
    output:
        work_dir+"/processed_info/{stem,[^_]+}_isIg.tsv"
    params:
        hmm = "data/hhms/immunoglobin/imuno.hmm"
    log: 
        out = work_dir+"/processed_info/{stem}_isIg.log"
    threads: 2
    shell:
        """
        hmmscan --tblout {output}  --noali   {params.hmm} {input} &> {log.out}
        """



rule determine_interacting_groups:
    input:
        ig_data = work_dir+"/processed_info/{stem}_isIg.tsv",
        fa = work_dir+"/processed_info/{stem,[^_]+}.fasta"
    output:
        work_dir+"/processed_info/{stem,[^_]+}_interactigGroups.tsv"
    params:
        predefined = "data/binding_partners.csv"
    threads: 2
    notebook:
        "notebooks/detect_interactors.r.ipynb"


# pdb file name structure
# {pdbname}_{frameid}_{part}.pdb
# {pdbname} - is a short name of PDB structure, cannot contain '_';
# {frameid} - a model id. In case of a static case keep it equal to '0', must be an integer;
# {part} - any of 'part1', 'part2', 'full'
rule split_complex:
    input: 
        pdb = work_dir+"/processed/{stem}.pdb",
        interacting_groups = work_dir+"/processed_info/{stem}_interactigGroups.tsv"
    output:
        part1 = work_dir+"/static/splits/{stem,[^_]+}_0_part1.pdb",
        part2 = work_dir+"/static/splits/{stem,[^_]+}_0_part2.pdb",
        full = work_dir+"/static/splits/{stem,[^_]+}_0_full.pdb",
    run:
        part1,part2 = get_interacting_chains(input.interacting_groups)
        part1sel=','.join(part1)
        part2sel=','.join(part2)
        print(part1sel,part2sel)
        shell(
            """
                pdb_selchain -{part1sel} {input.pdb}  > {output.part1} 
                pdb_selchain -{part2sel} {input.pdb} > {output.part2}
                cp {input.pdb}   {output.full} 
                
            """
        )
    
