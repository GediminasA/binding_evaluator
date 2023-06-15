### CONTAINER BUILDING RULES ###
rule singularity_container:
    input:
        "{container}.def"
    output:
        "{container}.sif"
    shell:
        "apptainer  build  {output} {input}"




def get_interacting_chains(precalculated_data):
    data = open(precalculated_data,"r").readlines()[0].replace("\n","").split("\t")
    data = [l.split(',') for l in data]
    data = [data[0],data[1]]
    data = sorted(data, key=lambda l: (len(l)),reverse=False)
    main_chain = data[0]
    out = data[1]
    out_str = " ".join(out)
    return(main_chain,list(out))
    
def get_interacting_chainsv1(pdbf,main_chain):
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
    print(f"Largest chain is {largest_chain_id}")
    is_main_chain = (main_chain in ids)
    chain_ids_str = " ".join(ids)
    if is_main_chain:
        print(f"The user indicated main chain {main_chain} is among the detected chains {chain_ids_str}")
    else:
        print(f"The user indicated main chain {main_chain} is not among the detected chains {chain_ids_str}")
    out=list()
    if not is_main_chain:
        main_chain = largest_chain_id
    out = set(ids)-set([main_chain])
    out_str = " ".join(out)
    print(f" Main chain is {main_chain}, other chains are: {out_str}")
    return(main_chain,list(out))
    