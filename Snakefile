import os
configfile: "configs/test.yaml"

pdb_stems = config["pdb_stems"].split()
work_dir = config["workdir"]




# include rules to prepare file downloadable rom file
include: "snakefiles/prepare_from_pdbdownloaded.smk"


rule initialfix:
    input:
        expand(work_dir+"/processed/{stem}.pdb",stem=pdb_stems)
