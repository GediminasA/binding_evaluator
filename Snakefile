import os
configfile: "configs/test.yaml"

pdb_stems = config["pdb_stems"].split()
work_dir = config["workdir"]




# include rules to prepare file downloadable rom file
include: "snakefiles/prepare_from_pdbdownloaded.smk"
# rosetta related evaluation
include: "snakefiles/rosetta_binding_evaluation.smk"

local_docking_runs_ids = list(range(1,config['rosetta']['local_docking_runs'],1))
rule rosetta_static_evaluation:
    input:
        #expand(work_dir+"/rosetta_on_static/{stem}_0001.pdb",stem=pdb_stems)
        #expand(work_dir+"/rosetta_on_static/localdockref_eval/{stem}_{id}.sc",stem=pdb_stems,id=local_docking_runs_ids)
        expand(work_dir+"/rosetta_on_static/localdockref/{stem}_{id}_prodigy.sc",stem=pdb_stems,id=local_docking_runs_ids)
        #expand(work_dir+"/scores/rosetta_static_{stem}_alldata.csv",stem=pdb_stems)

rule initialfix:
    input:
        expand(work_dir+"/processed/{stem}.pdb",stem=pdb_stems)
