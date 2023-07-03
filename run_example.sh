snakemake -c 6 t4t --use-conda --configfile configs/prodigy.yaml  --rerun-incomplete -k
snakemake -c 6 --use-conda --configfile configs/prodigy.yaml  --edit-notebook prodigy_run/initial_cleanup/1B6C_noambi.pdb
snakemake -c 2 --use-conda --configfile configs/prodigy.yaml    --use-conda --use-singularity -f get_freesasa
snakemake --profile local --configfile configs/antibody.yaml -f  mutants_target
snakemake --profile local --configfile configs/antibody.yaml -k collect_data_befor_modelling
snakemake --profile local --configfile configs/antibody.yaml -f antibody_run/mutants_structure_generation/EVOEF/structures/7LQV=A=E175K,F176R.pdb
snakemake --profile local --configfile configs/antibody.yaml -f  mutants_targets_templates
snakemake --profile local --configfile configs/antibody.yaml  -f  antibody_run/rezults/promod_models_results_main.csv
