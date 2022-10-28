import os
configfile: "configs/test.yaml"

pdb_stems = config["pdb_stems"].split()
work_dir = config["workdir"]



# include utilities 
include: "snakefiles/functions_and_utilities.smk"
# include rules to prepare file downloadable rom file
include: "snakefiles/prepare_from_pdbdownloaded.smk"
# prodigy evaluation on intial structure
include: "snakefiles/complex_evaluations_prodigy.smk"
# rosetta related evaluation
include: "snakefiles/rosetta_binding_evaluation.smk"
# solubility
include: "snakefiles/solubility_evaluations.smk"
# openmm forces
include: "snakefiles/openmm_forces.smk"


# main rules 

rule get_freesasa:
    input:
        expand(work_dir + "/static/splits/{stem}_0_{part}_freesasa.tsv",stem=pdb_stems,part=["part1","part2","full"])

rule get_splits:
    input:
        expand(work_dir + "/static/splits/{stem}_0_{part}.pdb",stem=pdb_stems,part=["part1","part2","full"])


rule prodigy_static_evaluation:
    input:
        expand(work_dir + "/static/{stem}_prodigy.tsv",stem=pdb_stems)

rule collect_data:
    input:
        expand(work_dir+"/processed_info/{stem}_interactigGroups.tsv", stem=pdb_stems)


rule calculate_openmm_ff1:
    input:
        expand(work_dir + "/static/splits/{stem}_0_{part}_ff_{forcefield}.tsv",stem=pdb_stems,part=["part1","part2","full"],forcefield=["amber99sbildn","amber10","amoeba2013","charmm36"])
