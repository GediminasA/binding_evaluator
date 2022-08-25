nice -n 19 snakemake --use-conda --conda-frontend mamba --configfile configs/test.yaml -c 80 --use-singularity -f  rosetta_static_evaluation
