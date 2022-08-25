import os
configfile: "configs/test.yaml"

pdb_stems = config["pdb_stems"].split()
work_dir = config["workdir"]




# include rules to prepare file downloadable rom file
include: "snakefiles/prepare_from_pdbdownloaded.smk"
# prodigy evaluation on intial structure
include: "snakefiles/complex_evaluations_prodigy.smk"
# rosetta related evaluation
include: "snakefiles/rosetta_binding_evaluation.smk"

rule prodigy_static_evaluation:
    input:
        expand(work_dir + "/prodigy_static/{stem}.sc",stem=pdb_stems)

rule rosetta_static_evaluation:
    input:
        expand(work_dir+"/scores/rosetta_static_{stem}_summary.csv",stem=pdb_stems)

rule initialfix:
    input:
        expand(work_dir+"/processed/{stem}.pdb",stem=pdb_stems)

rule collect_data:
    input:
        expand(work_dir+"/processed_info/{stem}_isIg.tsv", stem=pdb_stems)
