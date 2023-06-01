snakemake -c 6 t4t --use-conda --configfile configs/prodigy.yaml  --rerun-incomplete -k
snakemake -c 6 --use-conda --configfile configs/prodigy.yaml  --edit-notebook prodigy_run/initial_cleanup/1B6C_noambi.pdb
snakemake -c 2 --use-conda --configfile configs/prodigy.yaml    --use-conda --use-singularity -f get_freesasa
