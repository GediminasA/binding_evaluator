import os
configfile: "configs/test.yaml"


pdb_stems = config["pdb_stems"].split()
work_dir = config["workdir"]


# returns chains remaining  without a main chain
# if the expected main chain doesn't exits then the largest chain is taken as the main chain
def get_interacting_chains(pdbf,main_chain):
    from Bio import PDB as pdb
    ids = list()
    p = pdb.PDBParser()
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

    print("WW ",largest_chain_id)
    is_main_chain = (main_chain in ids)
    out=list()
    if not is_main_chain:
        main_chain = largest_chain_id
    out = set(ids)-set([main_chain])

    return(main_chain,list(out))






rule all:
    input: 
        expand(work_dir+"/pdbfixed/{stem}.pdb",stem=pdb_stems),
        expand(
            work_dir + "/pdbfixed_prodigy/{id}.sc",
            id=pdb_stems
        )

rule download_pdb:
    output:
        config["structures_folder"]+"/{stem}.pdb"
    params:
        stem = "{stem}"
    shell:
        """
            pdb_fetch {params.stem} > {output}
        """




rule copy_pdb_4_analysis:
    input:
        config["structures_folder"]+"/{stem}.pdb"
    output:
        work_dir + "/pristine/{stem}.pdb"
    run:
        out_full = os.path.abspath(output[0]) 
        in_full = os.path.abspath(input[0]) 
        shell(
        """
            cp   {in_full} {out_full}
        """
        )
    

rule remove_hetatoms:
    input:
        work_dir + "/pristine/{stem}.pdb"
    output:
        work_dir + "/delhetatm/{stem}.pdb"
    shell:
        """
            pdb_delhetatm {input} > {output}
        """


rule pdb_fix:
    input:
        work_dir + "/delhetatm/{stem}.pdb"
    output:
        work_dir + "/pdbfixed/{stem}.pdb"
    shell:
        """
        pdbfixer {input} --output {output} --add-residues
        """




rule run_prodigy_pdbfixed:
    input:
        complex_str = work_dir + "/pdbfixed/{id}.pdb"
    output:
        complex_prodigy = work_dir + "/pdbfixed_prodigy/{id}.sc"
    params:
        id= "{id}",
        main_chain=config["main_chain"]
    run:
        print(params.id)
        part1,part2 = get_interacting_chains(input[0],params.main_chain)
        part2 = ",".join(part2)
        shell(
            """
                set +e
                out=`prodigy  {input} --selection {part1} {part2} -q `
                exitcode=$?
                echo $exitcode
                if [ $exitcode -eq 1 ]
                then
                    echo NA  {part1} {part2} {params.id} > {output}
                else
                    echo $out | cut -f 2 -d " " | sed "s/$/ {part1} {part2} {params.id}/g" > {output} 
                fi
            """
        )
