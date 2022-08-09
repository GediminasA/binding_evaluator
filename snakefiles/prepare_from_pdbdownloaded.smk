import os



rule singularity_container:
    input:
        "container.def"
    output:
        "container.simg"
    shell:
        "singularity  build --fakeroot {output} {input}"


# returns chains remaining  without a main chain
# if the expected main chain doesn't exits then the largest chain is taken as the main chain
def get_interacting_chains(pdbf,main_chain):
    from Bio import PDB as pdb
    ids = list()
    p = pdb.PDBParser(QUIET=True)
    structure = p.get_structure("X", pdbf)
    sizes = []
    for model in structure:
        for chain in model:
            residues_n = len(list(chain.get_residues()))
            id = chain.get_id()
            ids.append(id)
            sizes.append((id,residues_n))
    sizes.sort(key=lambda tup: tup[1],reverse = True)
    largest_chain_id = sizes[0][0]

    is_main_chain = (main_chain in ids)
    out=list()
    if not is_main_chain:
        main_chain = largest_chain_id
    out = set(ids)-set([main_chain])

    return(main_chain,list(out))

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
        covid-lt/bin/pdb_align {input} 1>  {output} 2> {log}
        """ 

rule fix_with_promod:
    input:
        "{stem}.pdb",
        "container.simg"
    output:
        "{stem}_promod.pdb"
    log:
        "{stem}_promod.log"
    singularity: "container.simg"
    shell:
        """
        covid-lt/bin/promod-fix-pdb   {input[0]} 1>  {output} 2> {log}
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