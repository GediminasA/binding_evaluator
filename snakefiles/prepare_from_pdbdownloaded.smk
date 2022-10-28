import os






# returns chains remaining  without a main chain
# if the expected main chain doesn't exits then the largest chain is taken as the main chain

pdbproc_dir = work_dir + "/pdb_proc"

rule download_pdb:
    output:
        pdbproc_dir+"/pristine/{stem}.pdb"
        #config["structures_folder"]+"/{stem}.pdb"
    params:
        stem = "{stem}"
    container: None 
    shell:
        """
            pdb_fetch {params.stem} > {output}
        """


rule renumber_against_seqres:
    input:
        pdbproc_dir  + "/pristine/{stem}.pdb"
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
    shell:
        """
        PYTHONPATH=covid-lt covid-lt/bin/promod-fix-pdb {input[0]} 1> {output} 2> {log}
        """ 

rule pdb2pqr:
    input:
        "{stem}.pdb"
    output:
        pdb = "{stem}_pdb2pqr.pdb",
        pqr = "{stem}_pdb2pqr.pqr"
    log:
        "{stem}_pdb2pqr.log",
    shell:
        """
        pdb2pqr30 --drop-water  {input} {output.pqr}  --pdb-output {output.pdb} &> {log}
        """

rule pdb_fix:
    input:
        "{stem}.pdb"
    output:
        "{stem}_pdbfix.pdb"
    shell:
        """
        pdbfixer {input} --output {output} --add-residues
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
            work_dir+"/process_provided/{stem}_pdb2pqr_pdbfix.pdb"
        output:
            work_dir+"/processed/{stem,[^_]+}.pdb"
        shell:
            "cp {input} {output}"

rule copy_prepared:
        input:
            pdbproc_dir + "/process/{stem}_seqresMatched_promod_pdb2pqr_pdbfix.pdb"
        output:
            work_dir+"/processed/{stem,[^_]+}.pdb"
        shell:
            "cp {input} {output}"

ruleorder:  copy_alreadyprepared > copy_prepared

#extract fasta files and run hhblitz to identify antibodies

rule extract_seqs:
    input:
        work_dir+"/processed/{stem}.pdb"
    output:
        work_dir+"/processed_info/{stem,[^_]+}.fasta"
    shell:
        "pdb_tofasta -multi {input} > {output} "

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
    
